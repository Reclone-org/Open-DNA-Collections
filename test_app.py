#!/usr/bin/env python3
"""
Test script for the Open DNA Collections Streamlit app
This script tests the core functionality without running the full Streamlit interface
"""

import pandas as pd
import os
import glob
from pathlib import Path
from Bio import SeqIO
try:
    from Bio.SeqUtils import GC
except ImportError:
    # For newer versions of BioPython
    from Bio.SeqUtils import gc_fraction
    def GC(seq):
        return gc_fraction(seq) * 100
import sys

def test_data_loading():
    """Test data loading functionality."""
    print("🔍 Testing data loading...")
    
    base_path = "/Users/adaravena/Documents/GIT_PROJECTS/Open-DNA-Collections"
    
    # Test main CSV loading
    main_csv_path = Path(base_path) / "odc_plasmids.csv"
    if main_csv_path.exists():
        df = pd.read_csv(main_csv_path)
        print(f"✅ Main CSV loaded: {len(df)} records")
        print(f"   Columns: {list(df.columns)}")
        print(f"   Collections: {df['Collection'].nunique()}")
    else:
        print("❌ Main CSV not found")
        return False
    
    # Test collection-specific CSVs
    csv_files = [
        "Ecoli Nanobody Toolkit/ecoli_nb_tkit.csv",
        "Open Yeast Collection/Platemaps/OYC-v1_0.csv"
    ]
    
    for csv_file in csv_files:
        full_path = Path(base_path) / csv_file
        if full_path.exists():
            try:
                collection_df = pd.read_csv(full_path)
                print(f"✅ {csv_file}: {len(collection_df)} records")
            except Exception as e:
                print(f"⚠️  Error loading {csv_file}: {e}")
        else:
            print(f"⚠️  {csv_file} not found")
    
    # Test GenBank files
    genbank_pattern = str(Path(base_path) / "genbank" / "*.gb")
    gb_files = glob.glob(genbank_pattern)
    
    if gb_files:
        print(f"✅ Found {len(gb_files)} GenBank files")
        
        # Test loading a few GenBank files
        test_count = min(3, len(gb_files))
        for i, gb_file in enumerate(gb_files[:test_count]):
            try:
                with open(gb_file, 'r') as f:
                    record = SeqIO.read(f, "genbank")
                    gc_content = GC(record.seq)
                    print(f"   - {record.id}: {len(record.seq)} bp, GC: {gc_content:.1f}%")
            except Exception as e:
                print(f"   - Error reading {gb_file}: {e}")
    else:
        print("⚠️  No GenBank files found")
    
    return True

def test_search_functionality():
    """Test search functionality."""
    print("\n🔍 Testing search functionality...")
    
    base_path = "/Users/adaravena/Documents/GIT_PROJECTS/Open-DNA-Collections"
    main_csv_path = Path(base_path) / "odc_plasmids.csv"
    
    if not main_csv_path.exists():
        print("❌ Cannot test search - main CSV not available")
        return False
    
    df = pd.read_csv(main_csv_path)
    
    # Test various search scenarios
    test_queries = [
        "DNA polymerase",
        "BBF10K_003247",
        "ODC_0007",
        "Open Enzyme"
    ]
    
    for query in test_queries:
        mask = (
            df['Name'].str.contains(query, case=False, na=False) |
            df['BBF ID'].str.contains(query, case=False, na=False) |
            df['ODC ID'].str.contains(query, case=False, na=False) |
            df['Collection'].str.contains(query, case=False, na=False)
        )
        results = df[mask]
        print(f"✅ Query '{query}': {len(results)} results")
    
    return True

def test_analytics_data():
    """Test analytics data preparation."""
    print("\n📊 Testing analytics data...")
    
    base_path = "/Users/adaravena/Documents/GIT_PROJECTS/Open-DNA-Collections"
    main_csv_path = Path(base_path) / "odc_plasmids.csv"
    
    if not main_csv_path.exists():
        print("❌ Cannot test analytics - main CSV not available")
        return False
    
    df = pd.read_csv(main_csv_path)
    
    # Collection distribution
    collection_counts = df['Collection'].value_counts()
    print(f"✅ Collection distribution calculated:")
    for collection, count in collection_counts.head().items():
        print(f"   - {collection}: {count}")
    
    # ID completion
    bbf_complete = df['BBF ID'].notna().sum()
    odc_complete = df['ODC ID'].notna().sum()
    print(f"✅ ID completion: BBF: {bbf_complete}/{len(df)}, ODC: {odc_complete}/{len(df)}")
    
    return True

def test_imports():
    """Test if all required packages can be imported."""
    print("📦 Testing package imports...")
    
    required_packages = [
        'streamlit',
        'pandas', 
        'plotly',
        'Bio',
        'pathlib',
        'glob'
    ]
    
    for package in required_packages:
        try:
            if package == 'Bio':
                from Bio import SeqIO
                try:
                    from Bio.SeqUtils import GC
                except ImportError:
                    from Bio.SeqUtils import gc_fraction
            else:
                __import__(package)
            print(f"✅ {package}")
        except ImportError as e:
            print(f"❌ {package}: {e}")
            return False
    
    return True

def main():
    """Run all tests."""
    print("🧬 Open DNA Collections Streamlit App - Test Suite")
    print("=" * 60)
    
    all_tests_passed = True
    
    # Test imports
    if not test_imports():
        print("\n❌ Import tests failed. Please install required packages:")
        print("   pip install -r requirements.txt")
        all_tests_passed = False
    
    # Test data loading
    if not test_data_loading():
        all_tests_passed = False
    
    # Test search functionality
    if not test_search_functionality():
        all_tests_passed = False
    
    # Test analytics data
    if not test_analytics_data():
        all_tests_passed = False
    
    print("\n" + "=" * 60)
    if all_tests_passed:
        print("✅ All tests passed! The Streamlit app should work correctly.")
        print("\nTo run the app:")
        print("   streamlit run streamlit_app.py")
    else:
        print("❌ Some tests failed. Please check the issues above.")
    
    return all_tests_passed

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
