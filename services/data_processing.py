"""Data normalization and cache-building helpers for Open DNA Collections."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

MAIN_REQUIRED_COLUMNS = ["BBF ID", "ODC ID", "Name", "Collection"]

# Output columns used by the app and cache artifacts.
PLATEMAP_CANONICAL_COLUMNS = [
    "Platemap_Key",
    "Toolkit",
    "Platemap_Version",
    "Source_Path",
    "Well_Location",
    "Bacterial_Resistance",
    "Growth_Strain",
    "Growth_Conditions",
    "ODC ID",
    "BBF ID",
    "Name",
    "ODC_ID_NORM",
    "BBF_ID_NORM",
]

ENRICH_COLUMNS = [
    "Platemap_Key",
    "Toolkit",
    "Platemap_Version",
    "Well_Location",
    "Bacterial_Resistance",
    "Growth_Strain",
    "Growth_Conditions",
]


def _normalized_col_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _strip_bom(value: str) -> str:
    return str(value).replace("\ufeff", "")


def normalize_id(raw_value: object) -> Optional[str]:
    """Normalize ODC/BBF identifiers into stable lookup keys."""
    if raw_value is None:
        return None
    value = str(raw_value).strip()
    if not value or value.lower() == "nan":
        return None

    value = value.upper().replace(" ", "").replace("-", "_")
    value = re.sub(r"_+", "_", value)

    odc_match = re.match(r"^ODC_?(\d+)$", value)
    if odc_match:
        return f"ODC_{odc_match.group(1).zfill(4)}"

    bbf_match = re.match(r"^BBF10K_?(\d+)$", value)
    if bbf_match:
        return f"BBF10K_{bbf_match.group(1).zfill(6)}"

    return value


def _find_matching_column(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    keys = {_normalized_col_key(c): c for c in candidates}
    for col in df.columns:
        key = _normalized_col_key(_strip_bom(col))
        if key in keys:
            return col
    return None


def _extract_column(df: pd.DataFrame, candidates: Iterable[str]) -> pd.Series:
    col = _find_matching_column(df, candidates)
    if col is None:
        return pd.Series([None] * len(df), index=df.index, dtype="object")
    return df[col]


def normalize_main_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize main metadata CSV structure and add lookup key columns."""
    result = df.copy()
    result.columns = [_strip_bom(str(c)).strip() for c in result.columns]

    # Attempt to align column names if headers vary.
    rename_map: Dict[str, str] = {}
    alias_candidates = {
        "BBF ID": ["BBF ID", "bbf_id", "bbfid"],
        "ODC ID": ["ODC ID", "odc_id", "odcid"],
        "Name": ["Name", "part_name", "gene or insert name"],
        "Collection": ["Collection", "toolkit", "parent collection"],
    }
    for canonical, aliases in alias_candidates.items():
        col = _find_matching_column(result, aliases)
        if col and col != canonical:
            rename_map[col] = canonical

    if rename_map:
        result = result.rename(columns=rename_map)

    for col in MAIN_REQUIRED_COLUMNS:
        if col not in result.columns:
            result[col] = None

    # Keep only meaningful text normalization for lookup stability.
    result["ODC_ID_NORM"] = result["ODC ID"].apply(normalize_id)
    result["BBF_ID_NORM"] = result["BBF ID"].apply(normalize_id)
    result["Name"] = result["Name"].fillna("").astype(str).str.strip()
    result["Collection"] = result["Collection"].fillna("").astype(str).str.strip()
    return result


def normalize_platemap_dataframe(
    df: pd.DataFrame,
    toolkit: str,
    version: str,
    source_path: str,
) -> pd.DataFrame:
    """Normalize heterogeneous platemap CSV columns into a canonical shape."""
    normalized = pd.DataFrame(index=df.index)
    normalized["Platemap_Key"] = f"{toolkit}_{version}"
    normalized["Toolkit"] = toolkit
    normalized["Platemap_Version"] = version
    normalized["Source_Path"] = source_path

    normalized["Well_Location"] = _extract_column(
        df,
        ["Well Location", "well_location", "well address", "well_address"],
    )
    normalized["Bacterial_Resistance"] = _extract_column(
        df,
        ["Bacterial Resistance", "bacterial resistance", "well_media"],
    )
    normalized["Growth_Strain"] = _extract_column(
        df,
        ["Growth Strain", "Expression Strain", "growth_strain", "expression_strain"],
    )
    normalized["Growth_Conditions"] = _extract_column(
        df,
        [
            "Growth Conditions",
            "E. coli Growth Conditions",
            "growth_conditions",
            "ecoli_growth_conditions",
            "well_media",
        ],
    )
    normalized["ODC ID"] = _extract_column(df, ["ODC ID", "odc_id", "odc id"])
    normalized["BBF ID"] = _extract_column(
        df,
        ["BBF ID", "bbf_id", "bbf id", "part_gene_id", "part gene id"],
    )
    normalized["Name"] = _extract_column(
        df,
        ["Name", "Gene or Insert Name", "part_name", "product type"],
    )

    normalized["ODC_ID_NORM"] = normalized["ODC ID"].apply(normalize_id)
    normalized["BBF_ID_NORM"] = normalized["BBF ID"].apply(normalize_id)

    # Ensure stable column order for cache artifacts.
    for col in PLATEMAP_CANONICAL_COLUMNS:
        if col not in normalized.columns:
            normalized[col] = None

    return normalized[PLATEMAP_CANONICAL_COLUMNS]


def build_platemap_lookups(platemaps_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Build hash-join lookup tables keyed by normalized ODC and BBF IDs."""
    if platemaps_df.empty:
        empty = pd.DataFrame(columns=ENRICH_COLUMNS)
        return empty, empty

    odc_lookup = (
        platemaps_df[platemaps_df["ODC_ID_NORM"].notna()]
        .drop_duplicates(subset=["ODC_ID_NORM"], keep="first")
        .set_index("ODC_ID_NORM")[ENRICH_COLUMNS]
    )

    bbf_lookup = (
        platemaps_df[platemaps_df["BBF_ID_NORM"].notna()]
        .drop_duplicates(subset=["BBF_ID_NORM"], keep="first")
        .set_index("BBF_ID_NORM")[ENRICH_COLUMNS]
    )

    return odc_lookup, bbf_lookup


def discover_local_platemap_files(base_path: Path) -> List[Path]:
    """Return all toolkit platemap files from a repository checkout."""
    return sorted(base_path.glob("*/Platemaps/*.csv"))
