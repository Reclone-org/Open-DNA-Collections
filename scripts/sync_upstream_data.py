#!/usr/bin/env python3
"""Sync metadata from upstream Open-DNA-Collections and build cache artifacts.

This script is designed for CI scheduling (GitHub Actions) and local execution.
It writes cache artifacts to data/cache/ for fast Streamlit runtime loading.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import logging
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
from urllib.parse import quote

import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from services.data_processing import (
    build_platemap_lookups,
    discover_local_platemap_files,
    normalize_id,
    normalize_main_dataframe,
    normalize_platemap_dataframe,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("sync_upstream_data")

SCHEMA_VERSION = "1.0.0"
UPSTREAM_DEFAULT_REPO = "Reclone-org/Open-DNA-Collections"
UPSTREAM_DEFAULT_BRANCH = "main"


@dataclass
class SourceMetadata:
    source_repo: str
    source_branch: str
    source_commit: str


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_head_sha(base_dir: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(base_dir),
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except Exception:
        return "unknown"


def _normalize_csv_text(text: str) -> str:
    # Ensure UTF-8 BOM doesn't leak into header names.
    return text.replace("\ufeff", "")


def _is_platemap_path(path: str) -> bool:
    return bool(re.match(r"^[^/]+/Platemaps/.+\.csv$", path))


class GitHubUpstreamClient:
    def __init__(self, repo: str, branch: str, token: Optional[str] = None):
        self.repo = repo
        self.branch = branch
        self.session = requests.Session()
        headers = {"Accept": "application/vnd.github+json"}
        if token:
            headers["Authorization"] = f"Bearer {token}"
        self.session.headers.update(headers)

    def _get_json(self, url: str, params: Optional[Dict] = None) -> Dict:
        response = self.session.get(url, params=params, timeout=60)
        response.raise_for_status()
        return response.json()

    def get_branch_commit(self) -> str:
        branch_url = f"https://api.github.com/repos/{self.repo}/branches/{self.branch}"
        data = self._get_json(branch_url)
        return data["commit"]["sha"]

    def list_target_files(self, commit_sha: str) -> List[str]:
        tree_url = f"https://api.github.com/repos/{self.repo}/git/trees/{commit_sha}"
        tree_data = self._get_json(tree_url, params={"recursive": "1"})

        matched: List[str] = []
        for node in tree_data.get("tree", []):
            if node.get("type") != "blob":
                continue
            path = node.get("path", "")
            if path == "odc_plasmids.csv" or _is_platemap_path(path):
                matched.append(path)
        return sorted(matched)

    def fetch_text_file(self, path: str) -> str:
        encoded_path = quote(path, safe="/")
        raw_url = f"https://raw.githubusercontent.com/{self.repo}/{self.branch}/{encoded_path}"
        response = self.session.get(raw_url, timeout=60)
        response.raise_for_status()
        return _normalize_csv_text(response.text)


def load_data_from_upstream(
    repo: str,
    branch: str,
    token: Optional[str],
) -> Tuple[pd.DataFrame, pd.DataFrame, SourceMetadata]:
    client = GitHubUpstreamClient(repo=repo, branch=branch, token=token)
    commit_sha = client.get_branch_commit()
    files = client.list_target_files(commit_sha)

    logger.info("Found %s target CSV files in upstream %s@%s", len(files), repo, branch)

    main_df = pd.DataFrame()
    platemap_frames: List[pd.DataFrame] = []

    for path in files:
        text = client.fetch_text_file(path)
        df = pd.read_csv(io.StringIO(text))

        if path == "odc_plasmids.csv":
            main_df = normalize_main_dataframe(df)
            continue

        if _is_platemap_path(path):
            toolkit = path.split("/", 1)[0]
            version = Path(path).stem
            platemap_frames.append(
                normalize_platemap_dataframe(
                    df=df,
                    toolkit=toolkit,
                    version=version,
                    source_path=path,
                )
            )

    platemaps_df = pd.concat(platemap_frames, ignore_index=True) if platemap_frames else pd.DataFrame()

    source = SourceMetadata(
        source_repo=repo,
        source_branch=branch,
        source_commit=commit_sha,
    )
    return main_df, platemaps_df, source


def load_data_from_local(base_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame, SourceMetadata]:
    main_path = base_dir / "odc_plasmids.csv"
    if not main_path.exists():
        raise FileNotFoundError(f"Missing local main CSV: {main_path}")

    main_df = normalize_main_dataframe(pd.read_csv(main_path))

    platemap_frames: List[pd.DataFrame] = []
    for platemap_path in discover_local_platemap_files(base_dir):
        df = pd.read_csv(platemap_path)
        toolkit = platemap_path.parent.parent.name
        version = platemap_path.stem
        platemap_frames.append(
            normalize_platemap_dataframe(
                df=df,
                toolkit=toolkit,
                version=version,
                source_path=str(platemap_path.relative_to(base_dir)),
            )
        )

    platemaps_df = pd.concat(platemap_frames, ignore_index=True) if platemap_frames else pd.DataFrame()

    source = SourceMetadata(
        source_repo="local-checkout",
        source_branch="local",
        source_commit=_git_head_sha(base_dir),
    )
    return main_df, platemaps_df, source


def _parse_genbank_file(file_path: Path) -> Optional[Dict]:
    """Parse minimal GenBank fields without loading heavy dependencies."""
    record_id = None
    description = ""
    in_definition = False
    in_origin = False
    seq_parts: List[str] = []

    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("LOCUS"):
                    parts = line.split()
                    if len(parts) > 1:
                        record_id = parts[1]
                elif line.startswith("DEFINITION"):
                    description = line.replace("DEFINITION", "", 1).strip()
                    in_definition = True
                elif in_definition:
                    if line.startswith(" "):
                        description += " " + line.strip()
                    else:
                        in_definition = False

                if line.startswith("ORIGIN"):
                    in_origin = True
                    continue
                if in_origin:
                    if line.startswith("//"):
                        break
                    seq_parts.append("".join(ch for ch in line if ch.isalpha()))

        sequence = "".join(seq_parts).upper()
        if not sequence:
            return None

        gc_count = sequence.count("G") + sequence.count("C")
        gc_content = (gc_count / len(sequence) * 100.0) if sequence else 0.0

        part_id = normalize_id(file_path.stem) or normalize_id(record_id) or file_path.stem

        return {
            "part_id": part_id,
            "file_path": str(file_path),
            "length": len(sequence),
            "gc_content": gc_content,
            "record_id": record_id,
            "description": description.strip(),
            "sequence": sequence,
        }
    except Exception as exc:
        logger.warning("Failed parsing GenBank %s: %s", file_path, exc)
        return None


def build_genbank_assets(base_dir: Path, cache_dir: Path) -> Tuple[pd.DataFrame, Path]:
    patterns = [
        "genbank/*.gb",
        "*/Plasmids_Genbank/*.gb",
        "*/genbank_seq/*.gb",
    ]

    files: List[Path] = []
    for pattern in patterns:
        files.extend(base_dir.glob(pattern))
    files = sorted(set(files))

    rows = []
    blast_fasta_path = cache_dir / "blast_local_sequences.fasta"

    with blast_fasta_path.open("w", encoding="utf-8") as fasta_out:
        for file_path in files:
            parsed = _parse_genbank_file(file_path)
            if not parsed:
                continue

            relative_path = str(file_path.relative_to(base_dir))
            row = {k: v for k, v in parsed.items() if k != "sequence"}
            row["file_path"] = relative_path
            rows.append(row)

            fasta_out.write(f">{parsed['part_id']}\n{parsed['sequence']}\n")

    genbank_df = pd.DataFrame(rows)
    if not genbank_df.empty:
        genbank_df = genbank_df.drop_duplicates(subset=["part_id"], keep="first")

    return genbank_df, blast_fasta_path


def write_artifacts(
    base_dir: Path,
    cache_dir: Path,
    main_df: pd.DataFrame,
    platemaps_df: pd.DataFrame,
    source_meta: SourceMetadata,
) -> Dict:
    cache_dir.mkdir(parents=True, exist_ok=True)

    odc_lookup, bbf_lookup = build_platemap_lookups(platemaps_df)

    main_path = cache_dir / "main.parquet"
    platemap_path = cache_dir / "platemaps.parquet"
    odc_lookup_path = cache_dir / "platemap_lookup_odc.parquet"
    bbf_lookup_path = cache_dir / "platemap_lookup_bbf.parquet"
    genbank_index_path = cache_dir / "genbank_index.parquet"

    main_df.to_parquet(main_path, index=False)
    platemaps_df.to_parquet(platemap_path, index=False)
    odc_lookup.reset_index().to_parquet(odc_lookup_path, index=False)
    bbf_lookup.reset_index().to_parquet(bbf_lookup_path, index=False)

    genbank_df, blast_fasta_path = build_genbank_assets(base_dir=base_dir, cache_dir=cache_dir)
    genbank_df.to_parquet(genbank_index_path, index=False)

    tracked_files = [
        main_path,
        platemap_path,
        odc_lookup_path,
        bbf_lookup_path,
        genbank_index_path,
        blast_fasta_path,
    ]

    files_meta = []
    for file_path in tracked_files:
        files_meta.append(
            {
                "path": str(file_path.relative_to(base_dir)),
                "size_bytes": file_path.stat().st_size if file_path.exists() else 0,
                "sha256": _sha256_file(file_path) if file_path.exists() else None,
            }
        )

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "source_repo": source_meta.source_repo,
        "source_branch": source_meta.source_branch,
        "source_commit": source_meta.source_commit,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "files": files_meta,
        "counts": {
            "main_rows": int(len(main_df)),
            "platemap_rows": int(len(platemaps_df)),
            "platemap_lookup_odc": int(len(odc_lookup)),
            "platemap_lookup_bbf": int(len(bbf_lookup)),
            "genbank_rows": int(len(genbank_df)),
        },
    }

    manifest_path = cache_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    logger.info("Wrote cache artifacts to %s", cache_dir)
    logger.info("Counts: %s", manifest["counts"])

    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sync upstream ODC data and build cache artifacts")
    parser.add_argument("--base-dir", default=".", help="Repository root path")
    parser.add_argument(
        "--source",
        choices=["upstream", "local"],
        default="upstream",
        help="Data source mode",
    )
    parser.add_argument("--repo", default=UPSTREAM_DEFAULT_REPO, help="Upstream GitHub repo owner/name")
    parser.add_argument("--branch", default=UPSTREAM_DEFAULT_BRANCH, help="Upstream branch to sync")
    parser.add_argument(
        "--token-env",
        default="GITHUB_TOKEN",
        help="Environment variable with GitHub token for API calls",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()
    cache_dir = base_dir / "data" / "cache"

    if args.source == "upstream":
        token = os.getenv(args.token_env)
        main_df, platemaps_df, source_meta = load_data_from_upstream(
            repo=args.repo,
            branch=args.branch,
            token=token,
        )
    else:
        main_df, platemaps_df, source_meta = load_data_from_local(base_dir=base_dir)

    if main_df.empty:
        raise RuntimeError("Main dataset is empty after sync")

    write_artifacts(
        base_dir=base_dir,
        cache_dir=cache_dir,
        main_df=main_df,
        platemaps_df=platemaps_df,
        source_meta=source_meta,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
