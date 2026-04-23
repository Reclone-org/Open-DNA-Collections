#!/usr/bin/env python3
"""Open DNA Collections Interactive Database (Streamlit, cache-first architecture)."""

from __future__ import annotations

import base64
import json
import logging
import time
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import plotly.express as px
import streamlit as st

from services import BlastService, DNACollectionDataService

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


st.set_page_config(
    page_title="Open DNA Collections Database",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
<style>
    .main-header {
        font-size: 2.8rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 1.2rem;
    }
    .sequence-box {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        font-family: 'Courier New', monospace;
        font-size: 0.9rem;
        word-break: break-all;
        border: 1px solid #dee2e6;
        white-space: pre-wrap;
    }
    .freshness-card {
        background: #eef6fb;
        border: 1px solid #cde6f5;
        border-radius: 8px;
        padding: 12px;
        margin-bottom: 12px;
    }
</style>
""",
    unsafe_allow_html=True,
)


@st.cache_resource
def load_data_service() -> DNACollectionDataService:
    base_path = str(Path(__file__).parent)
    return DNACollectionDataService(base_path)


@st.cache_resource
def load_blast_service() -> BlastService:
    base_path = str(Path(__file__).parent)
    ncbi_email = None

    # Prefer Streamlit secrets, then environment variable.
    try:
        if "NCBI_EMAIL" in st.secrets:
            ncbi_email = st.secrets["NCBI_EMAIL"]
    except Exception:
        ncbi_email = None

    if not ncbi_email:
        ncbi_email = st.session_state.get("ncbi_email") or None

    return BlastService(base_path=base_path, ncbi_email=ncbi_email)


def create_download_link(df: pd.DataFrame, filename: str) -> str:
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    return f'<a href="data:file/csv;base64,{b64}" download="{filename}">Download CSV</a>'


def format_sequence(sequence: str, line_length: int = 80) -> str:
    lines = []
    for i in range(0, len(sequence), line_length):
        line = sequence[i : i + line_length]
        lines.append(f"{i + 1:>8}: {line}")
    return "\n".join(lines)


def show_data_freshness(service: DNACollectionDataService) -> None:
    freshness = service.get_data_freshness()
    diagnostics = service.get_diagnostics()

    st.markdown("### Data Freshness")
    st.markdown(
        f"""
<div class="freshness-card">
<strong>Source:</strong> {freshness.source_repo}<br/>
<strong>Branch:</strong> {freshness.source_branch}<br/>
<strong>Commit:</strong> {freshness.source_commit[:12]}<br/>
<strong>Generated:</strong> {freshness.generated_at}<br/>
<strong>Schema:</strong> {freshness.schema_version}
</div>
""",
        unsafe_allow_html=True,
    )

    with st.expander("Diagnostics"):
        st.json(diagnostics)



def show_home_page(service: DNACollectionDataService) -> None:
    st.markdown("## Welcome to the Open DNA Collections Database")
    st.markdown(
        """
This cache-first Streamlit app provides fast search across Reclone Open DNA Collections,
with upstream-synced platemap metadata and BLAST support.
"""
    )

    summary = service.get_collections_summary()

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        total_parts = sum(info["count"] for key, info in summary.items() if key != "GenBank Files")
        st.metric("Total Parts", total_parts)

    with col2:
        total_collections = len([k for k in summary.keys() if k != "GenBank Files"])
        st.metric("Collections", total_collections)

    with col3:
        total_genbank = summary.get("GenBank Files", {}).get("count", 0)
        st.metric("GenBank Files", total_genbank)

    with col4:
        avg_length = summary.get("GenBank Files", {}).get("avg_length", 0)
        st.metric("Avg. Sequence Length", f"{avg_length:.0f} bp")

    show_data_freshness(service)

    st.markdown("## Collections Overview")
    for collection, info in summary.items():
        if collection == "GenBank Files":
            continue
        with st.expander(f"{collection} ({info['count']} parts)"):
            c1, c2, c3 = st.columns(3)
            c1.metric("Total Parts", info["count"])
            c2.metric("With BBF ID", info["with_bbf_id"])
            c3.metric("With ODC ID", info["with_odc_id"])



def show_search_page(service: DNACollectionDataService) -> None:
    st.markdown("## 🔍 Search & Browse Collections")

    col1, col2 = st.columns([2, 1])
    with col1:
        search_query = st.text_input(
            "Search parts (Name, BBF ID, ODC ID, Collection):",
            placeholder="e.g., DNA polymerase, BBF10K_003247, ODC_0007",
        )

    with col2:
        collections = ["All Collections"]
        if not service.main_df.empty:
            collections.extend(sorted(service.main_df["Collection"].dropna().unique()))
        selected_collection = st.selectbox("Filter by Collection:", collections)

    with st.expander("🧬 Advanced Platemap Filters"):
        pmap_df = service.platemaps_df
        c1, c2, c3 = st.columns(3)

        with c1:
            resistance_types = ["All"]
            if not pmap_df.empty and "Bacterial_Resistance" in pmap_df.columns:
                resistance_types.extend(sorted(x for x in pmap_df["Bacterial_Resistance"].dropna().unique()))
            selected_resistance = st.selectbox("Bacterial Resistance:", list(dict.fromkeys(resistance_types)))

        with c2:
            growth_strains = ["All"]
            if not pmap_df.empty and "Growth_Strain" in pmap_df.columns:
                growth_strains.extend(sorted(x for x in pmap_df["Growth_Strain"].dropna().unique()))
            selected_strain = st.selectbox("Growth Strain:", list(dict.fromkeys(growth_strains)))

        with c3:
            well_pattern = st.text_input("Well Location Pattern:", placeholder="e.g., A1, A*, *1")

    if st.button("Search") or search_query:
        filters: Dict[str, str] = {}
        if selected_resistance != "All":
            filters["resistance"] = selected_resistance
        if selected_strain != "All":
            filters["strain"] = selected_strain
        if well_pattern:
            filters["well_pattern"] = well_pattern

        t0 = time.perf_counter()
        results = service.search_parts(
            query=search_query,
            collection=selected_collection,
            platemap_filters=filters,
        )
        elapsed_ms = (time.perf_counter() - t0) * 1000

        st.caption(f"Search latency: {elapsed_ms:.1f} ms")

        if results.empty:
            st.info("No results found. Try adjusting your search terms.")
            return

        st.markdown(f"### Found {len(results)} results")

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total Results", len(results))
        m2.metric("With Platemap Info", int(results["Well_Location"].notna().sum()))
        m3.metric("Resistance Types", int(results["Bacterial_Resistance"].nunique(dropna=True)))
        m4.metric("Growth Strains", int(results["Growth_Strain"].nunique(dropna=True)))

        d1, d2 = st.columns([1, 1])
        with d1:
            show_full_table = st.checkbox("Show all columns", value=False)
        with d2:
            items_per_page = st.selectbox("Items per page:", [10, 25, 50, 100], index=1)

        total_pages = (len(results) - 1) // items_per_page + 1
        page_num = st.selectbox("Page:", range(1, total_pages + 1), index=0)

        start_idx = (page_num - 1) * items_per_page
        page_results = results.iloc[start_idx : start_idx + items_per_page]

        if show_full_table:
            st.dataframe(page_results, use_container_width=True)
        else:
            display_cols = [
                "BBF ID",
                "ODC ID",
                "Name",
                "Collection",
                "Platemap_Version",
                "Well_Location",
                "Bacterial_Resistance",
                "Growth_Strain",
            ]
            available_cols = [col for col in display_cols if col in page_results.columns]
            display_results = page_results[available_cols].copy()
            for col in ["Well_Location", "Bacterial_Resistance", "Growth_Strain", "Platemap_Version"]:
                if col in display_results.columns:
                    display_results[col] = display_results[col].fillna("No info")
            st.dataframe(display_results, use_container_width=True)

        st.markdown(create_download_link(results, "search_results.csv"), unsafe_allow_html=True)

    st.markdown("## 📋 Available Platemaps")
    platemap_summary_df = service.get_platemap_summary()

    if platemap_summary_df.empty:
        st.info("No platemap cache data found.")
        return

    st.dataframe(platemap_summary_df, use_container_width=True)

    with st.expander("🔬 Detailed Platemap Viewer"):
        keys = platemap_summary_df["Platemap_Key"].tolist()
        selected_key = st.selectbox("Select Platemap:", keys)
        if selected_key:
            data = service.get_platemap(selected_key)
            st.dataframe(data, use_container_width=True)
            st.markdown(create_download_link(data, f"{selected_key}_platemap.csv"), unsafe_allow_html=True)



def show_analytics_page(service: DNACollectionDataService) -> None:
    st.markdown("## 📊 Analytics Dashboard")

    df = service.main_df
    if df.empty:
        st.warning("No main dataset available for analytics.")
        return

    st.markdown("### Collection Distribution")
    collection_counts = df["Collection"].value_counts()

    c1, c2 = st.columns(2)
    with c1:
        fig_pie = px.pie(
            values=collection_counts.values,
            names=collection_counts.index,
            title="Parts Distribution by Collection",
        )
        st.plotly_chart(fig_pie, use_container_width=True)

    with c2:
        fig_bar = px.bar(
            x=collection_counts.index,
            y=collection_counts.values,
            title="Parts Count by Collection",
            labels={"x": "Collection", "y": "Number of Parts"},
        )
        fig_bar.update_layout(xaxis_tickangle=45)
        st.plotly_chart(fig_bar, use_container_width=True)

    st.markdown("### ID Completion Analysis")
    id_stats = pd.DataFrame(
        {
            "BBF ID": [df["BBF ID"].notna().sum(), df["BBF ID"].isna().sum()],
            "ODC ID": [df["ODC ID"].notna().sum(), df["ODC ID"].isna().sum()],
        },
        index=["Present", "Missing"],
    )
    st.plotly_chart(px.bar(id_stats, title="ID Completion Status", barmode="group"), use_container_width=True)

    seq_stats = service.get_sequence_stats()
    if seq_stats.empty:
        st.info("No sequence index data available for sequence analytics.")
        return

    st.markdown("### Sequence Analysis")
    a1, a2 = st.columns(2)

    with a1:
        fig_length = px.histogram(
            seq_stats,
            x="length",
            title="Sequence Length Distribution",
            labels={"length": "Sequence Length (bp)", "count": "Count"},
            nbins=25,
        )
        st.plotly_chart(fig_length, use_container_width=True)

    with a2:
        fig_gc = px.histogram(
            seq_stats,
            x="gc_content",
            title="GC Content Distribution",
            labels={"gc_content": "GC Content (%)", "count": "Count"},
            nbins=25,
        )
        st.plotly_chart(fig_gc, use_container_width=True)

    st.markdown("### Sequence Statistics")
    s1, s2, s3, s4 = st.columns(4)
    s1.metric("Total Sequences", len(seq_stats))
    s2.metric("Avg Length", f"{seq_stats['length'].mean():.0f} bp")
    s3.metric("Min Length", f"{seq_stats['length'].min():.0f} bp")
    s4.metric("Max Length", f"{seq_stats['length'].max():.0f} bp")



def show_part_details_page(service: DNACollectionDataService) -> None:
    st.markdown("## 📋 Part Details")

    part_id = st.text_input(
        "Enter Part ID (BBF ID or ODC ID):",
        placeholder="e.g., BBF10K_003247 or ODC_0007",
    )

    if not part_id:
        return

    details = service.get_part_details(part_id)
    if not details:
        st.warning(f"No details found for part ID: {part_id}")
        return

    if "basic_info" in details:
        st.markdown("### Basic Information")
        info = details["basic_info"]
        c1, c2 = st.columns(2)
        c1.info(f"**BBF ID:** {info.get('BBF ID', 'N/A')}")
        c1.info(f"**ODC ID:** {info.get('ODC ID', 'N/A')}")
        c2.info(f"**Name:** {info.get('Name', 'N/A')}")
        c2.info(f"**Collection:** {info.get('Collection', 'N/A')}")

    if "genbank" not in details:
        st.info("No GenBank file found for this part.")
        return

    gb_info = details["genbank"]
    st.markdown("### Sequence Information")

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Length", f"{gb_info.get('length', 0)} bp")
    m2.metric("GC Content", f"{gb_info.get('gc_content', 0):.1f}%")
    m3.metric("Features", gb_info.get("features", 0))

    with m4:
        file_path = gb_info.get("file_path")
        if file_path and Path(file_path).exists():
            file_text = Path(file_path).read_text(encoding="utf-8", errors="ignore")
            st.download_button("Download GenBank", file_text, f"{part_id}.gb", "text/plain")

    if gb_info.get("description"):
        st.markdown("**Description:**")
        st.write(gb_info["description"])

    sequence = gb_info.get("sequence", "")
    if not sequence:
        return

    st.markdown("### DNA Sequence")
    o1, o2 = st.columns(2)
    with o1:
        show_formatted = st.checkbox("Show formatted sequence", value=True)
    with o2:
        line_length = st.slider("Line length", 40, 120, 80, step=10)

    if show_formatted:
        st.markdown(
            f'<div class="sequence-box">{format_sequence(sequence, line_length)}</div>',
            unsafe_allow_html=True,
        )
    else:
        st.code(sequence)

    st.download_button(
        "Download FASTA",
        f">{part_id}\n{sequence}",
        f"{part_id}.fasta",
        "text/plain",
    )



def show_data_management_page(service: DNACollectionDataService) -> None:
    st.markdown("## 📁 Data Management")

    diagnostics = service.get_diagnostics()

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Cache Overview:**")
        st.write(f"- Main rows: {diagnostics['main_rows']}")
        st.write(f"- Platemap rows: {diagnostics['platemap_rows']}")
        st.write(f"- ODC lookup keys: {diagnostics['platemap_odc_lookup_keys']}")
        st.write(f"- BBF lookup keys: {diagnostics['platemap_bbf_lookup_keys']}")

    with c2:
        st.markdown("**Sequence Index:**")
        st.write(f"- GenBank indexed rows: {diagnostics['genbank_index_rows']}")
        st.write(f"- Cache directory: `{diagnostics['cache_dir']}`")

    st.markdown("### Export Data")
    export_collection = st.selectbox(
        "Select dataset to export:",
        ["Main", "Platemaps", "GenBank Index"],
    )

    if export_collection == "Main":
        export_df = service.main_df
        filename = "open_dna_collections_main.csv"
    elif export_collection == "Platemaps":
        export_df = service.platemaps_df
        filename = "open_dna_collections_platemaps.csv"
    else:
        export_df = service.genbank_index_df
        filename = "open_dna_collections_genbank_index.csv"

    st.download_button(
        "Download CSV",
        export_df.to_csv(index=False),
        filename,
        "text/csv",
    )

    st.markdown("### Cache Manifest")
    if service.manifest:
        st.json(service.manifest)
    else:
        st.warning("No cache manifest found. Run scripts/sync_upstream_data.py to generate artifacts.")



def show_blast_page(service: DNACollectionDataService) -> None:
    st.markdown("## 🧪 BLAST Search")
    st.markdown(
        """
Run sequence similarity search against local Open DNA collection sequences.
Optional NCBI fallback can be used when local BLAST has no hits.
"""
    )

    if "ncbi_email" not in st.session_state:
        st.session_state["ncbi_email"] = ""

    with st.expander("NCBI Configuration (for remote fallback)"):
        st.session_state["ncbi_email"] = st.text_input(
            "NCBI Email",
            value=st.session_state["ncbi_email"],
            help="Required by NCBI for remote BLAST fallback.",
        )

    blast_service = BlastService(
        base_path=str(Path(__file__).parent),
        ncbi_email=st.session_state["ncbi_email"] or None,
    )

    c1, c2, c3 = st.columns(3)
    with c1:
        mode = st.selectbox(
            "Execution Mode",
            ["hybrid", "local", "ncbi"],
            index=0,
            help="hybrid = local first, fallback to NCBI when local has no hits",
        )
    with c2:
        target = st.selectbox("Target", ["auto", "local", "ncbi"], index=0)
    with c3:
        max_hits = st.slider("Max Hits", min_value=5, max_value=100, value=25, step=5)

    p1, p2 = st.columns(2)
    with p1:
        evalue = st.number_input("E-value threshold", min_value=0.0, value=1e-5, format="%.1e")
    with p2:
        poll_interval = st.number_input(
            "NCBI poll interval (sec)",
            min_value=60,
            max_value=300,
            value=60,
            step=30,
            help="NCBI recommends no more than 1 status poll per minute for a RID.",
        )

    query_sequence = st.text_area(
        "Query Sequence (DNA or Protein)",
        height=180,
        placeholder="Paste sequence here",
    )

    if st.button("Run BLAST"):
        if not query_sequence.strip():
            st.error("Please provide a sequence.")
            return

        params = {
            "max_hits": int(max_hits),
            "evalue": float(evalue),
            "poll_interval_sec": int(poll_interval),
            "timeout_sec": 300,
        }

        with st.spinner("Running BLAST query..."):
            start = time.perf_counter()
            result = blast_service.run_blast(
                mode=mode,
                sequence=query_sequence,
                target=target,
                params=params,
            )
            elapsed = time.perf_counter() - start

        st.markdown("### BLAST Result")
        st.caption(f"Elapsed: {elapsed:.2f} s")

        status = result.get("status", "unknown")
        source = result.get("source", "n/a")
        program = result.get("program", "n/a")
        query_type = result.get("query_type", "n/a")
        job_id = result.get("job_id", "n/a")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Status", status)
        c2.metric("Source", source)
        c3.metric("Program", program)
        c4.metric("Query Type", query_type)
        st.text(f"Job ID: {job_id}")

        if result.get("error"):
            st.error(result["error"])

        hits = result.get("hits", [])
        if hits:
            hits_df = pd.DataFrame(hits)
            st.dataframe(hits_df, use_container_width=True)
            st.download_button(
                "Download Hits (CSV)",
                hits_df.to_csv(index=False),
                f"blast_hits_{job_id}.csv",
                "text/csv",
            )
        else:
            st.info("No hits found.")

        raw_ref = result.get("raw_output_ref")
        if raw_ref:
            raw_path = Path(str(raw_ref))
            if raw_path.exists():
                st.download_button(
                    "Download Raw Output",
                    raw_path.read_text(encoding="utf-8", errors="ignore"),
                    raw_path.name,
                    "text/plain",
                )



def main() -> None:
    st.markdown('<h1 class="main-header">🧬 Open DNA Collections Database</h1>', unsafe_allow_html=True)
    st.markdown("### Streamlit cache-first browser for Reclone Open DNA Collections")

    try:
        service = load_data_service()
    except Exception as exc:
        st.error(f"Failed to load data service: {exc}")
        st.stop()

    st.sidebar.title("Navigation")
    if "current_page" not in st.session_state:
        st.session_state.current_page = "🏠 Home"

    pages = [
        "🏠 Home",
        "🔍 Search & Browse",
        "📊 Analytics",
        "📋 Part Details",
        "🧪 BLAST Search",
        "📁 Data Management",
    ]

    for page_name in pages:
        if st.sidebar.button(
            page_name,
            use_container_width=True,
            type="primary" if st.session_state.current_page == page_name else "secondary",
        ):
            st.session_state.current_page = page_name
            st.rerun()

    page = st.session_state.current_page
    if page == "🏠 Home":
        show_home_page(service)
    elif page == "🔍 Search & Browse":
        show_search_page(service)
    elif page == "📊 Analytics":
        show_analytics_page(service)
    elif page == "📋 Part Details":
        show_part_details_page(service)
    elif page == "🧪 BLAST Search":
        show_blast_page(service)
    elif page == "📁 Data Management":
        show_data_management_page(service)


if __name__ == "__main__":
    main()
