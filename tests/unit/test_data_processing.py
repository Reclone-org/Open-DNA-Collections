import pandas as pd

from services.data_processing import (
    build_platemap_lookups,
    normalize_id,
    normalize_main_dataframe,
    normalize_platemap_dataframe,
)


def test_normalize_id_variants():
    assert normalize_id("odc7") == "ODC_0007"
    assert normalize_id("ODC-0337") == "ODC_0337"
    assert normalize_id("bbf10k-42") == "BBF10K_000042"
    assert normalize_id(None) is None
    assert normalize_id("  ") is None


def test_normalize_main_dataframe_adds_lookup_columns():
    df = pd.DataFrame(
        {
            "\ufeffBBF ID": ["BBF10K_000589"],
            "ODC ID": ["ODC_0001"],
            "Name": ["pOpen"],
            "Collection": ["Open Plasmids"],
        }
    )
    normalized = normalize_main_dataframe(df)

    assert "BBF ID" in normalized.columns
    assert "ODC_ID_NORM" in normalized.columns
    assert normalized.loc[0, "ODC_ID_NORM"] == "ODC_0001"
    assert normalized.loc[0, "BBF_ID_NORM"] == "BBF10K_000589"


def test_normalize_platemap_dataframe_handles_alias_columns():
    raw = pd.DataFrame(
        {
            "well_address": ["A1"],
            "part_name": ["Example Part"],
            "part_gene_id": ["BBF10K_000589"],
            "well_media": ["ampicillin"],
        }
    )

    normalized = normalize_platemap_dataframe(
        df=raw,
        toolkit="Open Enzyme Collection",
        version="OEC-v0_1-1",
        source_path="Open Enzyme Collection/Platemaps/OEC-v0_1-1.csv",
    )

    assert normalized.loc[0, "Well_Location"] == "A1"
    assert normalized.loc[0, "Name"] == "Example Part"
    assert normalized.loc[0, "BBF_ID_NORM"] == "BBF10K_000589"
    assert normalized.loc[0, "Bacterial_Resistance"] == "ampicillin"


def test_build_platemap_lookups_prefers_first_match():
    platemaps = pd.DataFrame(
        {
            "Platemap_Key": ["K1", "K2"],
            "Toolkit": ["T1", "T2"],
            "Platemap_Version": ["v1", "v2"],
            "Well_Location": ["A1", "B2"],
            "Bacterial_Resistance": ["amp", "kan"],
            "Growth_Strain": ["DH5a", "BL21"],
            "Growth_Conditions": ["LB", "LB"],
            "ODC ID": ["ODC_0007", "ODC_0007"],
            "BBF ID": ["BBF10K_003247", "BBF10K_003247"],
            "Name": ["N1", "N2"],
            "ODC_ID_NORM": ["ODC_0007", "ODC_0007"],
            "BBF_ID_NORM": ["BBF10K_003247", "BBF10K_003247"],
        }
    )

    odc_lookup, bbf_lookup = build_platemap_lookups(platemaps)

    assert odc_lookup.loc["ODC_0007", "Well_Location"] == "A1"
    assert bbf_lookup.loc["BBF10K_003247", "Well_Location"] == "A1"
