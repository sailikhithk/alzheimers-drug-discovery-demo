import pandas as pd
import numpy as np
import os

def create_aligned_dataset():
    print("Loading datasets...")
    
    # Load Descriptors (X)
    df_desc = pd.read_csv('descriptors_output.csv')
    print(f"Descriptors shape (raw): {df_desc.shape}")
    
    # Deduplicate descriptors
    # Since we verified that duplicates have 0 variance, we can safely keep the first occurrence
    df_desc = df_desc.drop_duplicates(subset=['Name'])
    print(f"Descriptors shape (deduplicated): {df_desc.shape}")
    
    # Load Bioactivity Data (Y source)
    df_bio = pd.read_csv('Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv')
    print(f"Bioactivity shape: {df_bio.shape}")
    
    # Check for IDs
    if 'Name' not in df_desc.columns:
        print("Error: 'Name' column missing in descriptors_output.csv")
        return
    if 'Molecule ChEMBL ID' not in df_bio.columns:
        print("Error: 'Molecule ChEMBL ID' column missing in bioactivity file")
        return

    # Merge on ID
    # We ideally want to preserve the bioactivity rows (1319) and attach descriptors to them
    print("Merging on ID...")
    df_merged = pd.merge(
        df_bio, # Left join on bioactivity to keep its structure/order if preferred, or inner
        df_desc, 
        left_on='Molecule ChEMBL ID', 
        right_on='Name', 
        how='inner'
    )
    
    print(f"Merged shape: {df_merged.shape}")
    
    if len(df_merged) != len(df_bio):
        print(f"Warning: Merged size ({len(df_merged)}) differs from source bioactivity size ({len(df_bio)}).")
    else:
        print("Success: Merged size matches original bioactivity data size.")
    
    # Calculate pIC50
    print("Calculating pIC50...")
    
    def calculate_pic50(input_df):
        pIC50 = []
        for i in input_df['Standard Value']:
            molar = i * (10**-9)  # Convert nM to M
            pIC50.append(-np.log10(molar))

        input_df['pIC50'] = pIC50
        return input_df

    df_merged = calculate_pic50(df_merged)
    
    # Verify pIC50 range
    print("pIC50 Summary:")
    print(df_merged['pIC50'].describe())
    
    # Save the aligned dataset
    output_filename = 'aligned_dataset.csv'
    df_merged.to_csv(output_filename, index=False)
    print(f"Saved aligned dataset to {output_filename}")

if __name__ == "__main__":
    create_aligned_dataset()
