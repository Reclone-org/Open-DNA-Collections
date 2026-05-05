import textwrap

from scripts.sync_upstream_data import load_data_from_upstream


MAIN_CSV = textwrap.dedent(
    """\
BBF ID,ODC ID,Name,Collection
BBF10K_003247,ODC_0007,9N7,Open Enzyme Collection
"""
)

PLATEMAP_CSV = textwrap.dedent(
    """\
Well Location,ODC ID,BBF ID,Bacterial Resistance,Growth Strain,Growth Conditions,Name
A1,ODC_0007,BBF10K_003247,ampicillin,DH5a,LB,9N7
"""
)


def test_load_data_from_upstream_with_mocked_client(monkeypatch):
    monkeypatch.setattr(
        "scripts.sync_upstream_data.GitHubUpstreamClient.get_branch_commit",
        lambda self: "abc123",
    )
    monkeypatch.setattr(
        "scripts.sync_upstream_data.GitHubUpstreamClient.list_target_files",
        lambda self, commit_sha: [
            "odc_plasmids.csv",
            "Open Enzyme Collection/Platemaps/OEC-v1_1.csv",
        ],
    )

    def _fake_fetch(self, path):
        if path == "odc_plasmids.csv":
            return MAIN_CSV
        return PLATEMAP_CSV

    monkeypatch.setattr(
        "scripts.sync_upstream_data.GitHubUpstreamClient.fetch_text_file",
        _fake_fetch,
    )

    main_df, platemaps_df, source = load_data_from_upstream(
        repo="Reclone-org/Open-DNA-Collections",
        branch="main",
        token=None,
    )

    assert len(main_df) == 1
    assert len(platemaps_df) == 1
    assert main_df.loc[0, "ODC_ID_NORM"] == "ODC_0007"
    assert platemaps_df.loc[0, "ODC_ID_NORM"] == "ODC_0007"
    assert source.source_commit == "abc123"
