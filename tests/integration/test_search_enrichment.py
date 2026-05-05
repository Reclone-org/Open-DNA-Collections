from pathlib import Path

import pandas as pd

from services.cache_service import DNACollectionDataService


def test_search_enrichment_works_with_local_fallback(tmp_path):
    # Build minimal local repo-like structure without cache manifest.
    main_df = pd.DataFrame(
        {
            "BBF ID": ["BBF10K_003247"],
            "ODC ID": ["ODC_0007"],
            "Name": ["9°N-7 DNA polymerase"],
            "Collection": ["Open Enzyme Collection"],
        }
    )
    main_df.to_csv(tmp_path / "odc_plasmids.csv", index=False)

    platemap_dir = tmp_path / "Open Enzyme Collection" / "Platemaps"
    platemap_dir.mkdir(parents=True)
    platemap_df = pd.DataFrame(
        {
            "Well Location": ["A1"],
            "ODC ID": ["ODC_0007"],
            "BBF ID": ["BBF10K_003247"],
            "Bacterial Resistance": ["ampicillin"],
            "Growth Strain": ["DH5a"],
            "Growth Conditions": ["LB"],
            "Name": ["9°N-7 DNA polymerase"],
        }
    )
    platemap_df.to_csv(platemap_dir / "OEC-v1_1.csv", index=False)

    service = DNACollectionDataService(base_path=str(tmp_path))
    result = service.search_parts(query="ODC_0007")

    assert len(result) == 1
    assert result.iloc[0]["Well_Location"] == "A1"
    assert result.iloc[0]["Bacterial_Resistance"] == "ampicillin"
    assert result.iloc[0]["Growth_Strain"] == "DH5a"
