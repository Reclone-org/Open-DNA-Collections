from pathlib import Path

from scripts.sync_upstream_data import _is_platemap_path, _parse_genbank_file


def test_is_platemap_path():
    assert _is_platemap_path("Open Enzyme Collection/Platemaps/OEC-v1_1.csv")
    assert not _is_platemap_path("Open Enzyme Collection/README.md")


def test_parse_genbank_file_extracts_sequence(tmp_path):
    gb = tmp_path / "ODC_0001.gb"
    gb.write_text(
        """
LOCUS       ODC_0001               12 bp    DNA     circular
DEFINITION  Example construct.
ORIGIN
        1 atgcgcatttaa
//
""".strip()
        + "\n",
        encoding="utf-8",
    )

    parsed = _parse_genbank_file(gb)
    assert parsed is not None
    assert parsed["part_id"] == "ODC_0001"
    assert parsed["length"] == 12
    assert parsed["sequence"] == "ATGCGCATTTAA"
