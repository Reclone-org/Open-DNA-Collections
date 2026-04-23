"""BLAST service layer with local-first and optional NCBI fallback."""

from __future__ import annotations

import csv
import logging
import re
import shutil
import subprocess
import tempfile
import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional

import requests

from .data_processing import normalize_id

logger = logging.getLogger(__name__)


class BlastService:
    """Run sequence BLAST queries using local BLAST+ and optional NCBI API fallback."""

    LOCAL_OUTFMT = "6 qseqid sseqid pident length evalue bitscore qcovs sstart send"

    def __init__(
        self,
        base_path: str,
        ncbi_email: Optional[str] = None,
        ncbi_tool: str = "open_dna_collections_streamlit",
    ):
        self.base_path = Path(base_path)
        self.cache_dir = self.base_path / "data" / "cache"
        self.blast_dir = self.cache_dir / "blast"
        self.blast_dir.mkdir(parents=True, exist_ok=True)

        self.local_fasta_path = self.cache_dir / "blast_local_sequences.fasta"
        self.local_db_prefix = self.blast_dir / "odc_local_nucl"

        self.ncbi_email = ncbi_email
        self.ncbi_tool = ncbi_tool

    @staticmethod
    def _clean_sequence(sequence: str) -> str:
        return re.sub(r"\s+", "", (sequence or "")).upper()

    @classmethod
    def detect_query_type(cls, sequence: str) -> str:
        seq = cls._clean_sequence(sequence)
        if not seq:
            raise ValueError("Sequence is empty")

        dna_pattern = re.compile(r"^[ACGTUN]+$")
        protein_pattern = re.compile(r"^[ABCDEFGHIKLMNPQRSTVWXYZ\*]+$")

        if dna_pattern.match(seq):
            return "dna"
        if protein_pattern.match(seq):
            return "protein"
        raise ValueError("Sequence contains invalid characters for DNA or protein query")

    def _blast_binary_available(self, binary_name: str) -> bool:
        return shutil.which(binary_name) is not None

    def local_blast_available(self) -> bool:
        required = ["makeblastdb", "blastn", "tblastn"]
        return all(self._blast_binary_available(name) for name in required)

    def _list_genbank_files(self) -> List[Path]:
        patterns = [
            "genbank/*.gb",
            "*/Plasmids_Genbank/*.gb",
            "*/genbank_seq/*.gb",
        ]

        files: List[Path] = []
        for pattern in patterns:
            files.extend(self.base_path.glob(pattern))
        return sorted(set(files))

    def _ensure_local_fasta(self) -> Path:
        if self.local_fasta_path.exists() and self.local_fasta_path.stat().st_size > 0:
            return self.local_fasta_path

        # Build a local FASTA from all repository GenBank files.
        try:
            from Bio import SeqIO
        except Exception as exc:
            raise RuntimeError("Biopython is required to build local BLAST FASTA") from exc

        genbank_files = self._list_genbank_files()
        if not genbank_files:
            raise RuntimeError("No GenBank files found for local BLAST database")

        with self.local_fasta_path.open("w", encoding="utf-8") as fasta_out:
            for gb_file in genbank_files:
                try:
                    with gb_file.open("r", encoding="utf-8", errors="ignore") as handle:
                        record = SeqIO.read(handle, "genbank")
                    sequence = str(record.seq).upper().strip()
                    if not sequence:
                        continue

                    part_id = normalize_id(gb_file.stem) or normalize_id(record.id) or gb_file.stem
                    fasta_out.write(f">{part_id}\n{sequence}\n")
                except Exception as exc:
                    logger.warning("Skipping GenBank file during FASTA build (%s): %s", gb_file, exc)

        if self.local_fasta_path.stat().st_size <= 0:
            raise RuntimeError("Failed to build local BLAST FASTA from GenBank files")
        return self.local_fasta_path

    def _ensure_local_db(self) -> Path:
        if not self.local_blast_available():
            raise RuntimeError(
                "Local BLAST binaries not available. Install ncbi-blast+ (e.g., via packages.txt on Streamlit)."
            )

        fasta = self._ensure_local_fasta()
        expected_files = [
            self.local_db_prefix.with_suffix(".nhr"),
            self.local_db_prefix.with_suffix(".nin"),
            self.local_db_prefix.with_suffix(".nsq"),
        ]

        needs_build = not all(path.exists() for path in expected_files)
        if not needs_build:
            db_mtime = min(path.stat().st_mtime for path in expected_files)
            needs_build = db_mtime < fasta.stat().st_mtime

        if needs_build:
            cmd = [
                "makeblastdb",
                "-in",
                str(fasta),
                "-dbtype",
                "nucl",
                "-out",
                str(self.local_db_prefix),
                "-title",
                "Open DNA Collections Local Nucleotide DB",
            ]
            logger.info("Building local BLAST DB: %s", " ".join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                raise RuntimeError(
                    f"makeblastdb failed ({result.returncode}): {result.stderr.strip() or result.stdout.strip()}"
                )

        return self.local_db_prefix

    def _parse_local_hits(self, tsv_path: Path) -> List[Dict]:
        hits: List[Dict] = []
        if not tsv_path.exists():
            return hits

        with tsv_path.open("r", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if len(row) < 9:
                    continue
                hits.append(
                    {
                        "subject_id": row[1],
                        "identity": float(row[2]),
                        "alignment_length": int(float(row[3])),
                        "evalue": float(row[4]),
                        "bitscore": float(row[5]),
                        "qcov": float(row[6]) if row[6] else None,
                        "sstart": int(float(row[7])),
                        "send": int(float(row[8])),
                    }
                )
        return hits

    def _run_local_blast(self, sequence: str, query_type: str, params: Optional[Dict] = None) -> Dict:
        params = params or {}
        db_prefix = self._ensure_local_db()

        max_hits = int(params.get("max_hits", 25))
        evalue = float(params.get("evalue", 1e-5))
        job_id = f"local-{uuid.uuid4().hex[:12]}"

        with tempfile.TemporaryDirectory(prefix="odc_blast_", dir=str(self.blast_dir)) as temp_dir:
            temp_dir_path = Path(temp_dir)
            query_path = temp_dir_path / "query.fasta"
            out_path = temp_dir_path / "results.tsv"

            cleaned = self._clean_sequence(sequence)
            query_path.write_text(f">query\n{cleaned}\n", encoding="utf-8")

            program = "blastn" if query_type == "dna" else "tblastn"
            cmd = [
                program,
                "-query",
                str(query_path),
                "-db",
                str(db_prefix),
                "-out",
                str(out_path),
                "-outfmt",
                self.LOCAL_OUTFMT,
                "-max_target_seqs",
                str(max_hits),
                "-evalue",
                str(evalue),
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                return {
                    "job_id": job_id,
                    "status": "failed",
                    "source": "local",
                    "program": program,
                    "query_type": query_type,
                    "hits": [],
                    "raw_output_ref": None,
                    "error": result.stderr.strip() or result.stdout.strip() or "Local BLAST execution failed",
                }

            hits = self._parse_local_hits(out_path)
            return {
                "job_id": job_id,
                "status": "completed",
                "source": "local",
                "program": program,
                "query_type": query_type,
                "hits": hits,
                "raw_output_ref": str(out_path),
            }

    def _parse_ncbi_submit(self, response_text: str) -> Dict[str, Optional[str]]:
        rid_match = re.search(r"RID\s*=\s*([A-Z0-9-]+)", response_text)
        rtoe_match = re.search(r"RTOE\s*=\s*(\d+)", response_text)
        return {
            "rid": rid_match.group(1) if rid_match else None,
            "rtoe": rtoe_match.group(1) if rtoe_match else None,
        }

    def _parse_ncbi_tabular_hits(self, text: str) -> List[Dict]:
        hits: List[Dict] = []
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            row = line.split("\t")
            if len(row) < 12:
                continue
            # BLAST tabular default fields:
            # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            hits.append(
                {
                    "subject_id": row[1],
                    "identity": float(row[2]),
                    "alignment_length": int(float(row[3])),
                    "evalue": float(row[10]),
                    "bitscore": float(row[11]),
                    "qcov": None,
                    "sstart": int(float(row[8])),
                    "send": int(float(row[9])),
                }
            )
        return hits

    def _run_ncbi_blast(self, sequence: str, query_type: str, params: Optional[Dict] = None) -> Dict:
        params = params or {}
        if not self.ncbi_email:
            return {
                "job_id": f"ncbi-{uuid.uuid4().hex[:12]}",
                "status": "failed",
                "source": "ncbi",
                "program": "blastn" if query_type == "dna" else "blastp",
                "query_type": query_type,
                "hits": [],
                "raw_output_ref": None,
                "error": "NCBI fallback requires an email (set NCBI_EMAIL).",
            }

        max_hits = int(params.get("max_hits", 25))
        evalue = float(params.get("evalue", 1e-5))
        poll_interval_sec = max(int(params.get("poll_interval_sec", 60)), 60)
        timeout_sec = int(params.get("timeout_sec", 240))

        program = "blastn" if query_type == "dna" else "blastp"
        database = params.get("ncbi_database") or ("core_nt" if query_type == "dna" else "swissprot")

        blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        submit_payload = {
            "CMD": "Put",
            "PROGRAM": program,
            "DATABASE": database,
            "QUERY": self._clean_sequence(sequence),
            "HITLIST_SIZE": str(max_hits),
            "EXPECT": str(evalue),
            "tool": self.ncbi_tool,
            "email": self.ncbi_email,
        }

        submit_response = requests.post(blast_url, data=submit_payload, timeout=60)
        submit_response.raise_for_status()
        submit_meta = self._parse_ncbi_submit(submit_response.text)
        rid = submit_meta.get("rid")

        if not rid:
            return {
                "job_id": f"ncbi-{uuid.uuid4().hex[:12]}",
                "status": "failed",
                "source": "ncbi",
                "program": program,
                "query_type": query_type,
                "hits": [],
                "raw_output_ref": None,
                "error": "NCBI did not return a RID for submitted BLAST job.",
            }

        initial_wait = int(submit_meta.get("rtoe") or 0)
        if initial_wait > 0:
            time.sleep(min(initial_wait, poll_interval_sec))

        start = time.time()
        status_payload = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_OBJECT": "SearchInfo",
            "tool": self.ncbi_tool,
            "email": self.ncbi_email,
        }

        while True:
            status_response = requests.get(blast_url, params=status_payload, timeout=60)
            status_response.raise_for_status()
            status_text = status_response.text

            if "Status=WAITING" in status_text:
                if time.time() - start > timeout_sec:
                    return {
                        "job_id": rid,
                        "status": "timeout",
                        "source": "ncbi",
                        "program": program,
                        "query_type": query_type,
                        "hits": [],
                        "raw_output_ref": rid,
                        "error": "Timed out waiting for NCBI BLAST result.",
                    }
                time.sleep(poll_interval_sec)
                continue

            if "Status=FAILED" in status_text:
                return {
                    "job_id": rid,
                    "status": "failed",
                    "source": "ncbi",
                    "program": program,
                    "query_type": query_type,
                    "hits": [],
                    "raw_output_ref": rid,
                    "error": "NCBI BLAST job failed.",
                }

            if "Status=UNKNOWN" in status_text:
                return {
                    "job_id": rid,
                    "status": "failed",
                    "source": "ncbi",
                    "program": program,
                    "query_type": query_type,
                    "hits": [],
                    "raw_output_ref": rid,
                    "error": "NCBI BLAST job RID was not recognized.",
                }

            if "Status=READY" in status_text:
                has_hits = "ThereAreHits=yes" in status_text
                if not has_hits:
                    return {
                        "job_id": rid,
                        "status": "completed",
                        "source": "ncbi",
                        "program": program,
                        "query_type": query_type,
                        "hits": [],
                        "raw_output_ref": rid,
                    }

                result_payload = {
                    "CMD": "Get",
                    "RID": rid,
                    "FORMAT_TYPE": "Text",
                    "ALIGNMENT_VIEW": "Tabular",
                    "DESCRIPTIONS": str(max_hits),
                    "ALIGNMENTS": str(max_hits),
                    "tool": self.ncbi_tool,
                    "email": self.ncbi_email,
                }
                result_response = requests.get(blast_url, params=result_payload, timeout=60)
                result_response.raise_for_status()
                hits = self._parse_ncbi_tabular_hits(result_response.text)

                raw_output_path = self.blast_dir / f"ncbi_{rid}.txt"
                raw_output_path.write_text(result_response.text, encoding="utf-8")

                return {
                    "job_id": rid,
                    "status": "completed",
                    "source": "ncbi",
                    "program": program,
                    "query_type": query_type,
                    "hits": hits,
                    "raw_output_ref": str(raw_output_path),
                }

            # If status output is unexpected, avoid hot looping.
            time.sleep(poll_interval_sec)

    def run_blast(
        self,
        mode: str,
        sequence: str,
        target: str = "auto",
        params: Optional[Dict] = None,
    ) -> Dict:
        """Contract: run_blast(mode, sequence, target, params) -> job/result schema."""
        params = params or {}
        query_type = self.detect_query_type(sequence)

        normalized_mode = (mode or "hybrid").strip().lower()
        normalized_target = (target or "auto").strip().lower()

        if normalized_mode not in {"local", "hybrid", "ncbi"}:
            normalized_mode = "hybrid"

        allow_local = normalized_mode in {"local", "hybrid"} and normalized_target in {"auto", "local"}
        allow_ncbi = normalized_mode in {"hybrid", "ncbi"} and normalized_target in {"auto", "ncbi"}

        local_result = None
        if allow_local:
            try:
                local_result = self._run_local_blast(sequence=sequence, query_type=query_type, params=params)
            except Exception as exc:
                local_result = {
                    "job_id": f"local-{uuid.uuid4().hex[:12]}",
                    "status": "failed",
                    "source": "local",
                    "program": "blastn" if query_type == "dna" else "tblastn",
                    "query_type": query_type,
                    "hits": [],
                    "raw_output_ref": None,
                    "error": str(exc),
                }

            # Local-only mode always returns local result.
            if normalized_mode == "local":
                return local_result

            # Hybrid mode returns local if useful.
            if local_result.get("status") == "completed" and local_result.get("hits"):
                return local_result

        if allow_ncbi:
            try:
                return self._run_ncbi_blast(sequence=sequence, query_type=query_type, params=params)
            except Exception as exc:
                return {
                    "job_id": f"ncbi-{uuid.uuid4().hex[:12]}",
                    "status": "failed",
                    "source": "ncbi",
                    "program": "blastn" if query_type == "dna" else "blastp",
                    "query_type": query_type,
                    "hits": [],
                    "raw_output_ref": None,
                    "error": str(exc),
                }

        if local_result is not None:
            return local_result

        return {
            "job_id": f"blast-{uuid.uuid4().hex[:12]}",
            "status": "failed",
            "source": "local",
            "program": "blastn" if query_type == "dna" else "tblastn",
            "query_type": query_type,
            "hits": [],
            "raw_output_ref": None,
            "error": "No BLAST target available for requested mode/target combination.",
        }
