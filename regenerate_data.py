#!/usr/bin/env python3
"""
Script to regenerate the dataset with improved NaN handling.
This will fetch fresh data and ensure we have a full dataset without NaN abstracts.
"""

import subprocess
import sys
import os

def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"\n{'='*50}")
    print(f"Running: {description}")
    print(f"Command: {cmd}")
    print('='*50)
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {description}:")
        print(f"Exit code: {e.returncode}")
        print(f"Error output: {e.stderr}")
        return False

def main():
    print("Regenerating dataset with improved NaN handling...")
    
    # Step 1: Fetch fresh data with supplementation
    if not run_command("python data_ingestion/fetch_pubmed.py", "Fetching fresh PubMed data"):
        print("Failed to fetch data. Exiting.")
        sys.exit(1)
    
    # Step 2: Preprocess abstracts
    if not run_command("python preprocessing/preprocess_abstracts.py", "Preprocessing abstracts"):
        print("Failed to preprocess abstracts. Exiting.")
        sys.exit(1)
    
    # Step 3: Extract categories
    if not run_command("python preprocessing/extract_categories.py", "Extracting categories"):
        print("Failed to extract categories. Exiting.")
        sys.exit(1)
    
    # Step 4: Verify results
    print("\n" + "="*50)
    print("Verifying results...")
    print("="*50)
    
    try:
        import pandas as pd
        df = pd.read_csv('preprocessing/cleaned_abstracts.csv')
        print(f"Total articles: {len(df)}")
        print(f"Articles with NaN abstracts: {df['Abstract'].isna().sum()}")
        print(f"Articles with valid abstracts: {len(df) - df['Abstract'].isna().sum()}")
        
        if df['Abstract'].isna().sum() == 0:
            print("✅ SUCCESS: All articles have valid abstracts!")
        else:
            print(f"⚠️  WARNING: {df['Abstract'].isna().sum()} articles still have NaN abstracts")
            
    except Exception as e:
        print(f"Error verifying results: {e}")
    
    print("\nDataset regeneration complete!")

if __name__ == '__main__':
    main() 