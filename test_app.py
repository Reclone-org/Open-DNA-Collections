#!/usr/bin/env python3
"""Lightweight smoke tests for Open DNA Collections app services."""

from __future__ import annotations

import sys
from pathlib import Path

from services import BlastService, DNACollectionDataService


def smoke_data_service(base_path: Path) -> bool:
    print("🔍 Testing data service...")
    service = DNACollectionDataService(str(base_path))

    diagnostics = service.get_diagnostics()
    print(f"✅ Main rows: {diagnostics['main_rows']}")
    print(f"✅ Platemap rows: {diagnostics['platemap_rows']}")

    if diagnostics["main_rows"] <= 0:
        print("❌ No main dataset rows available")
        return False

    results = service.search_parts("ODC_0007")
    print(f"✅ Search ODC_0007 -> {len(results)} rows")

    return True


def smoke_blast_contract(base_path: Path) -> bool:
    print("\n🧪 Testing BLAST contract...")
    blast = BlastService(str(base_path))

    try:
        query_type = blast.detect_query_type("ATGCATGCATGC")
        print(f"✅ Query type detection: {query_type}")
    except Exception as exc:
        print(f"❌ Query type detection failed: {exc}")
        return False

    result = blast.run_blast(mode="local", sequence="ATGCATGCATGC", target="local")
    required_keys = {"job_id", "status", "source", "program", "query_type", "hits", "raw_output_ref"}
    missing = required_keys - set(result.keys())

    if missing:
        print(f"❌ BLAST result missing keys: {sorted(missing)}")
        return False

    print(f"✅ BLAST contract response status: {result['status']}")
    return True


def main() -> int:
    base_path = Path(__file__).resolve().parent
    print("🧬 Open DNA Collections - Smoke Test")
    print("=" * 50)

    ok_data = smoke_data_service(base_path)
    ok_blast = smoke_blast_contract(base_path)

    print("\n" + "=" * 50)
    if ok_data and ok_blast:
        print("✅ Smoke tests passed")
        print("\nRun app:")
        print("  streamlit run streamlit_app.py")
        return 0

    print("❌ Smoke tests failed")
    return 1


if __name__ == "__main__":
    sys.exit(main())
