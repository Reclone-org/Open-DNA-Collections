from services.blast_service import BlastService


def test_detect_query_type_dna_and_protein(tmp_path):
    svc = BlastService(base_path=str(tmp_path))
    assert svc.detect_query_type("ATGCATGC") == "dna"
    assert svc.detect_query_type("MTEYKLVVVG") == "protein"


def test_run_blast_local_failure_returns_contract(tmp_path):
    svc = BlastService(base_path=str(tmp_path))

    result = svc.run_blast(
        mode="local",
        sequence="ATGCATGCATGC",
        target="local",
        params={"max_hits": 5},
    )

    assert "job_id" in result
    assert result["source"] == "local"
    assert result["query_type"] == "dna"
    assert "hits" in result
    assert isinstance(result["hits"], list)
