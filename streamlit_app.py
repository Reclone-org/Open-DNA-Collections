#!/usr/bin/env python3
"""
Open DNA Collections Interactive Database
A Streamlit application for browsing and exploring the Open DNA Collections.

Author: GitHub Copilot
Date: August 4, 2025
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import glob
from pathlib import Path
import re
from Bio import SeqIO
try:
    from Bio.SeqUtils import GC
except ImportError:
    # For newer versions of BioPython
    from Bio.SeqUtils import gc_fraction
    def GC(seq):
        return gc_fraction(seq) * 100
import io
import base64
from typing import Dict, List, Optional, Tuple
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Page configuration
st.set_page_config(
    page_title="Open DNA Collections Database",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 2rem;
    }
    .collection-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
    .metric-card {
        background-color: #ffffff;
        padding: 1rem;
        border-radius: 0.5rem;
        border: 1px solid #e0e0e0;
        text-align: center;
    }
    .sequence-box {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        font-family: 'Courier New', monospace;
        font-size: 0.9rem;
        word-break: break-all;
        border: 1px solid #dee2e6;
    }
</style>
""", unsafe_allow_html=True)

class DNACollectionDatabase:
    """Main class for handling DNA collection data."""
    
    def __init__(self, base_path: str):
        self.base_path = Path(base_path)
        self.data = {}
        self.genbank_files = {}
        self._load_all_data()
    
    def _load_all_data(self):
        """Load all CSV data and GenBank files."""
        try:
            # Load main CSV
            main_csv_path = self.base_path / "odc_plasmids.csv"
            if main_csv_path.exists():
                self.data['main'] = pd.read_csv(main_csv_path)
                logger.info(f"Loaded main CSV with {len(self.data['main'])} records")
            
            # Load collection-specific CSVs
            self._load_collection_csvs()
            
            # Load GenBank files
            self._load_genbank_files()
            
        except Exception as e:
            logger.error(f"Error loading data: {e}")
            st.error(f"Error loading data: {e}")
    
    def _load_collection_csvs(self):
        """Load collection-specific CSV files."""
        csv_patterns = [
            "Ecoli Nanobody Toolkit/ecoli_nb_tkit.csv",
            "Open Yeast Collection/Platemaps/OYC-v1_0.csv",
            "product-csvs/*.csv"
        ]
        
        for pattern in csv_patterns:
            csv_files = glob.glob(str(self.base_path / pattern))
            for csv_file in csv_files:
                try:
                    collection_name = Path(csv_file).stem
                    self.data[collection_name] = pd.read_csv(csv_file)
                    logger.info(f"Loaded {collection_name} with {len(self.data[collection_name])} records")
                except Exception as e:
                    logger.warning(f"Could not load {csv_file}: {e}")
    
    def _load_genbank_files(self):
        """Load GenBank files."""
        genbank_patterns = [
            "genbank/*.gb",
            "*/Plasmids_Genbank/*.gb",
            "*/genbank_seq/*.gb"
        ]
        
        for pattern in genbank_patterns:
            gb_files = glob.glob(str(self.base_path / pattern))
            for gb_file in gb_files:
                try:
                    with open(gb_file, 'r') as f:
                        record = SeqIO.read(f, "genbank")
                        self.genbank_files[record.id] = {
                            'record': record,
                            'file_path': gb_file,
                            'sequence': str(record.seq),
                            'length': len(record.seq),
                            'gc_content': GC(record.seq),
                            'features': len(record.features),
                            'description': record.description
                        }
                except Exception as e:
                    logger.warning(f"Could not load GenBank file {gb_file}: {e}")
        
        logger.info(f"Loaded {len(self.genbank_files)} GenBank files")
    
    def get_collections_summary(self) -> Dict:
        """Get summary statistics for all collections."""
        summary = {}
        
        if 'main' in self.data:
            df = self.data['main']
            collections = df['Collection'].value_counts()
            
            for collection, count in collections.items():
                summary[collection] = {
                    'count': count,
                    'with_bbf_id': df[df['Collection'] == collection]['BBF ID'].notna().sum(),
                    'with_odc_id': df[df['Collection'] == collection]['ODC ID'].notna().sum()
                }
        
        summary['GenBank Files'] = {
            'count': len(self.genbank_files),
            'total_sequences': len(self.genbank_files),
            'avg_length': sum(info['length'] for info in self.genbank_files.values()) / len(self.genbank_files) if self.genbank_files else 0
        }
        
        return summary
    
    def search_parts(self, query: str, collection: Optional[str] = None) -> pd.DataFrame:
        """Search for parts across all collections."""
        results = []
        
        if 'main' in self.data:
            df = self.data['main'].copy()
            
            if collection and collection != "All Collections":
                df = df[df['Collection'] == collection]
            
            if query:
                mask = (
                    df['Name'].str.contains(query, case=False, na=False) |
                    df['BBF ID'].str.contains(query, case=False, na=False) |
                    df['ODC ID'].str.contains(query, case=False, na=False) |
                    df['Collection'].str.contains(query, case=False, na=False)
                )
                df = df[mask]
            
            results.append(df)
        
        return pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    
    def get_part_details(self, part_id: str) -> Dict:
        """Get detailed information about a specific part."""
        details = {}
        
        # Search in main data
        if 'main' in self.data:
            match = self.data['main'][
                (self.data['main']['BBF ID'] == part_id) | 
                (self.data['main']['ODC ID'] == part_id)
            ]
            if not match.empty:
                details['basic_info'] = match.iloc[0].to_dict()
        
        # Search in GenBank files
        if part_id in self.genbank_files:
            details['genbank'] = self.genbank_files[part_id]
        
        return details

# Initialize the database
@st.cache_resource
def load_database():
    """Load the database with caching."""
    # Use relative path for deployment
    base_path = str(Path(__file__).parent)
    return DNACollectionDatabase(base_path)

def create_download_link(df: pd.DataFrame, filename: str) -> str:
    """Create a download link for a DataFrame."""
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">Download CSV</a>'
    return href

def format_sequence(sequence: str, line_length: int = 80) -> str:
    """Format a DNA sequence for display."""
    formatted = []
    for i in range(0, len(sequence), line_length):
        line = sequence[i:i+line_length]
        # Add position numbers
        pos = f"{i+1:>8}: "
        formatted.append(pos + line)
    return "\n".join(formatted)

def main():
    """Main Streamlit application."""
    
    # Header
    st.markdown('<h1 class="main-header">🧬 Open DNA Collections Database</h1>', unsafe_allow_html=True)
    st.markdown("### Interactive browser for the Reclone Open DNA Collections")
    
    # Load database
    try:
        db = load_database()
    except Exception as e:
        st.error(f"Failed to load database: {e}")
        st.stop()
    
    # Sidebar navigation
    st.sidebar.title("Navigation")
    
    # Initialize session state for current page
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "🏠 Home"
    
    # Navigation buttons instead of selectbox
    pages = ["🏠 Home", "🔍 Search & Browse", "📊 Analytics", "📋 Part Details", "📁 Data Management"]
    
    for page_name in pages:
        if st.sidebar.button(page_name, use_container_width=True, 
                           type="primary" if st.session_state.current_page == page_name else "secondary"):
            st.session_state.current_page = page_name
            st.rerun()
    
    # Display current page
    page = st.session_state.current_page
    
    # Main content based on selected page
    if page == "🏠 Home":
        show_home_page(db)
    elif page == "🔍 Search & Browse":
        show_search_page(db)
    elif page == "📊 Analytics":
        show_analytics_page(db)
    elif page == "📋 Part Details":
        show_part_details_page(db)
    elif page == "📁 Data Management":
        show_data_management_page(db)

def show_home_page(db: DNACollectionDatabase):
    """Display the home page with overview."""
    
    st.markdown("## Welcome to the Open DNA Collections Database")
    st.markdown("""
    This interactive database provides access to the Reclone Open DNA Collections, 
    a global collaboration for equitable access to biotechnology. Explore thousands of 
    DNA parts, plasmids, and biological tools across multiple collections.
    """)
    
    # Quick stats
    summary = db.get_collections_summary()
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        total_parts = sum(info['count'] for info in summary.values() if 'count' in info)
        st.metric("Total Parts", total_parts)
    
    with col2:
        total_collections = len([k for k in summary.keys() if k != 'GenBank Files'])
        st.metric("Collections", total_collections)
    
    with col3:
        total_genbank = summary.get('GenBank Files', {}).get('count', 0)
        st.metric("GenBank Files", total_genbank)
    
    with col4:
        avg_length = summary.get('GenBank Files', {}).get('avg_length', 0)
        st.metric("Avg. Sequence Length", f"{avg_length:.0f} bp")
    
    # Collections overview
    st.markdown("## Collections Overview")
    
    for collection, info in summary.items():
        if collection != 'GenBank Files':
            with st.expander(f"{collection} ({info['count']} parts)"):
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total Parts", info['count'])
                with col2:
                    st.metric("With BBF ID", info['with_bbf_id'])
                with col3:
                    st.metric("With ODC ID", info['with_odc_id'])
    
    # Quick Actions
    st.markdown("## Quick Actions")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("🔍 Browse All Parts", use_container_width=True):
            st.session_state.current_page = "🔍 Search & Browse"
            st.rerun()
    
    with col2:
        if st.button("📊 View Analytics", use_container_width=True):
            st.session_state.current_page = "📊 Analytics"
            st.rerun()
    
    with col3:
        if st.button("📁 Manage Data", use_container_width=True):
            st.session_state.current_page = "📁 Data Management"
            st.rerun()

def show_search_page(db: DNACollectionDatabase):
    """Display the search and browse page."""
    
    st.markdown("## 🔍 Search & Browse Collections")
    
    # Search controls
    col1, col2 = st.columns([2, 1])
    
    with col1:
        search_query = st.text_input(
            "Search parts (Name, BBF ID, ODC ID, Collection):",
            placeholder="e.g., DNA polymerase, BBF10K_003247, ODC_0007"
        )
    
    with col2:
        collections = ["All Collections"]
        if 'main' in db.data:
            collections.extend(sorted(db.data['main']['Collection'].dropna().unique()))
        
        selected_collection = st.selectbox("Filter by Collection:", collections)
    
    # Search results
    if st.button("Search") or search_query:
        results = db.search_parts(search_query, selected_collection)
        
        if not results.empty:
            st.markdown(f"### Found {len(results)} results")
            
            # Display options
            col1, col2 = st.columns([1, 1])
            with col1:
                show_full_table = st.checkbox("Show all columns", value=False)
            with col2:
                items_per_page = st.selectbox("Items per page:", [10, 25, 50, 100], index=1)
            
            # Pagination
            total_pages = (len(results) - 1) // items_per_page + 1
            page_num = st.selectbox("Page:", range(1, total_pages + 1), index=0)
            
            start_idx = (page_num - 1) * items_per_page
            end_idx = start_idx + items_per_page
            page_results = results.iloc[start_idx:end_idx]
            
            # Display results
            if show_full_table:
                st.dataframe(page_results, use_container_width=True)
            else:
                # Display key columns only
                display_cols = ['BBF ID', 'ODC ID', 'Name', 'Collection']
                available_cols = [col for col in display_cols if col in page_results.columns]
                st.dataframe(page_results[available_cols], use_container_width=True)
            
            # Download option
            st.markdown(create_download_link(results, "search_results.csv"), unsafe_allow_html=True)
            
        else:
            st.info("No results found. Try adjusting your search terms.")
    
    # Browse by collection
    st.markdown("## Browse by Collection")
    
    if 'main' in db.data:
        collections_data = db.data['main']['Collection'].value_counts()
        
        for collection, count in collections_data.items():
            with st.expander(f"{collection} ({count} parts)"):
                collection_data = db.data['main'][db.data['main']['Collection'] == collection]
                st.dataframe(collection_data[['BBF ID', 'ODC ID', 'Name']], use_container_width=True)

def show_analytics_page(db: DNACollectionDatabase):
    """Display analytics and visualizations."""
    
    st.markdown("## 📊 Analytics Dashboard")
    
    if 'main' not in db.data:
        st.warning("No main dataset available for analytics.")
        return
    
    df = db.data['main']
    
    # Collection distribution
    st.markdown("### Collection Distribution")
    collection_counts = df['Collection'].value_counts()
    
    col1, col2 = st.columns(2)
    
    with col1:
        fig_pie = px.pie(
            values=collection_counts.values,
            names=collection_counts.index,
            title="Parts Distribution by Collection"
        )
        st.plotly_chart(fig_pie, use_container_width=True)
    
    with col2:
        fig_bar = px.bar(
            x=collection_counts.index,
            y=collection_counts.values,
            title="Parts Count by Collection",
            labels={'x': 'Collection', 'y': 'Number of Parts'}
        )
        fig_bar.update_layout(xaxis_tickangle=45)
        st.plotly_chart(fig_bar, use_container_width=True)
    
    # ID completion analysis
    st.markdown("### ID Completion Analysis")
    
    id_stats = pd.DataFrame({
        'BBF ID': [df['BBF ID'].notna().sum(), df['BBF ID'].isna().sum()],
        'ODC ID': [df['ODC ID'].notna().sum(), df['ODC ID'].isna().sum()]
    }, index=['Present', 'Missing'])
    
    fig_id = px.bar(
        id_stats,
        title="ID Completion Status",
        barmode='group'
    )
    st.plotly_chart(fig_id, use_container_width=True)
    
    # GenBank file analysis
    if db.genbank_files:
        st.markdown("### Sequence Analysis")
        
        # Sequence length distribution
        lengths = [info['length'] for info in db.genbank_files.values()]
        gc_contents = [info['gc_content'] for info in db.genbank_files.values()]
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig_length = px.histogram(
                x=lengths,
                title="Sequence Length Distribution",
                labels={'x': 'Sequence Length (bp)', 'y': 'Count'},
                nbins=20
            )
            st.plotly_chart(fig_length, use_container_width=True)
        
        with col2:
            fig_gc = px.histogram(
                x=gc_contents,
                title="GC Content Distribution",
                labels={'x': 'GC Content (%)', 'y': 'Count'},
                nbins=20
            )
            st.plotly_chart(fig_gc, use_container_width=True)
        
        # Summary statistics
        st.markdown("### Sequence Statistics")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Sequences", len(lengths))
        with col2:
            st.metric("Avg Length", f"{sum(lengths)/len(lengths):.0f} bp")
        with col3:
            st.metric("Min Length", f"{min(lengths)} bp")
        with col4:
            st.metric("Max Length", f"{max(lengths)} bp")

def show_part_details_page(db: DNACollectionDatabase):
    """Display detailed part information."""
    
    st.markdown("## 📋 Part Details")
    
    # Part selection
    part_id = st.text_input(
        "Enter Part ID (BBF ID or ODC ID):",
        placeholder="e.g., BBF10K_003247 or ODC_0007"
    )
    
    if part_id:
        details = db.get_part_details(part_id)
        
        if details:
            # Basic information
            if 'basic_info' in details:
                st.markdown("### Basic Information")
                info = details['basic_info']
                
                col1, col2 = st.columns(2)
                with col1:
                    st.info(f"**BBF ID:** {info.get('BBF ID', 'N/A')}")
                    st.info(f"**ODC ID:** {info.get('ODC ID', 'N/A')}")
                with col2:
                    st.info(f"**Name:** {info.get('Name', 'N/A')}")
                    st.info(f"**Collection:** {info.get('Collection', 'N/A')}")
            
            # GenBank information
            if 'genbank' in details:
                st.markdown("### Sequence Information")
                gb_info = details['genbank']
                
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Length", f"{gb_info['length']} bp")
                with col2:
                    st.metric("GC Content", f"{gb_info['gc_content']:.1f}%")
                with col3:
                    st.metric("Features", gb_info['features'])
                with col4:
                    if st.button("Download GenBank"):
                        # Create download for GenBank file
                        with open(gb_info['file_path'], 'r') as f:
                            st.download_button(
                                "Download",
                                f.read(),
                                f"{part_id}.gb",
                                "text/plain"
                            )
                
                # Description
                if gb_info.get('description'):
                    st.markdown("**Description:**")
                    st.write(gb_info['description'])
                
                # Sequence display
                st.markdown("### DNA Sequence")
                
                # Sequence options
                col1, col2 = st.columns(2)
                with col1:
                    show_formatted = st.checkbox("Show formatted sequence", value=True)
                with col2:
                    line_length = st.slider("Line length", 40, 120, 80, step=10)
                
                if show_formatted:
                    formatted_seq = format_sequence(gb_info['sequence'], line_length)
                    st.markdown(f'<div class="sequence-box">{formatted_seq}</div>', 
                              unsafe_allow_html=True)
                else:
                    st.code(gb_info['sequence'], language=None)
                
                # Download sequence
                st.download_button(
                    "Download FASTA",
                    f">{part_id}\n{gb_info['sequence']}",
                    f"{part_id}.fasta",
                    "text/plain"
                )
            
        else:
            st.warning(f"No details found for part ID: {part_id}")

def show_data_management_page(db: DNACollectionDatabase):
    """Display data management tools."""
    
    st.markdown("## 📁 Data Management")
    
    # Data overview
    st.markdown("### Database Overview")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Loaded Datasets:**")
        for name, df in db.data.items():
            st.write(f"- {name}: {len(df)} records")
    
    with col2:
        st.markdown("**GenBank Files:**")
        st.write(f"- Total files: {len(db.genbank_files)}")
        st.write(f"- Total sequences: {len(db.genbank_files)}")
    
    # Export options
    st.markdown("### Export Data")
    
    export_collection = st.selectbox(
        "Select dataset to export:",
        ["All Data"] + list(db.data.keys())
    )
    
    if st.button("Export CSV"):
        if export_collection == "All Data":
            # Combine all data
            all_data = []
            for name, df in db.data.items():
                df_copy = df.copy()
                df_copy['source_dataset'] = name
                all_data.append(df_copy)
            
            if all_data:
                combined_df = pd.concat(all_data, ignore_index=True, sort=False)
                st.download_button(
                    "Download Combined CSV",
                    combined_df.to_csv(index=False),
                    "open_dna_collections_all.csv",
                    "text/csv"
                )
        else:
            if export_collection in db.data:
                df = db.data[export_collection]
                st.download_button(
                    f"Download {export_collection} CSV",
                    df.to_csv(index=False),
                    f"open_dna_collections_{export_collection}.csv",
                    "text/csv"
                )
    
    # Data validation
    st.markdown("### Data Validation")
    
    if st.button("Run Validation Checks"):
        validation_results = []
        
        if 'main' in db.data:
            df = db.data['main']
            
            # Check for missing data
            missing_bbf = df['BBF ID'].isna().sum()
            missing_odc = df['ODC ID'].isna().sum()
            missing_names = df['Name'].isna().sum()
            
            validation_results.extend([
                f"Missing BBF IDs: {missing_bbf}",
                f"Missing ODC IDs: {missing_odc}",
                f"Missing Names: {missing_names}"
            ])
            
            # Check for duplicates
            dup_bbf = df['BBF ID'].duplicated().sum()
            dup_odc = df['ODC ID'].duplicated().sum()
            
            validation_results.extend([
                f"Duplicate BBF IDs: {dup_bbf}",
                f"Duplicate ODC IDs: {dup_odc}"
            ])
        
        # Display results
        st.markdown("**Validation Results:**")
        for result in validation_results:
            st.write(f"- {result}")
    
    # File system info
    st.markdown("### File System Information")
    
    if st.button("Show File Locations"):
        st.markdown("**CSV Files:**")
        for name in db.data.keys():
            st.write(f"- {name}")
        
        st.markdown("**GenBank Files:**")
        gb_dirs = set()
        for info in db.genbank_files.values():
            gb_dirs.add(str(Path(info['file_path']).parent))
        
        for gb_dir in sorted(gb_dirs):
            file_count = sum(1 for info in db.genbank_files.values() 
                           if str(Path(info['file_path']).parent) == gb_dir)
            st.write(f"- {gb_dir}: {file_count} files")

if __name__ == "__main__":
    main()
