import pandas as pd
import os

def check_alignment():
    print("Loading files...")
    
    # Load descriptors
    if not os.path.exists('descriptors_output.csv'):
        print("Error: descriptors_output.csv not found")
        return
    df_desc = pd.read_csv('descriptors_output.csv')
    print(f"Descriptors shape: {df_desc.shape}")
    
    # Load source data (checking the curated and labelled file which contains IDs)
    if not os.path.exists('Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv'):
        print("Error: Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv not found")
        return
    df_source = pd.read_csv('Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv')
    print(f"Source data shape: {df_source.shape}")
    
    # Check lengths
    if len(df_desc) != len(df_source):
        print(f"WARNING: Length mismatch! Descriptors: {len(df_desc)}, Source: {len(df_source)}")
    else:
        print("Lengths match.")

    # Check ID alignment
    # 'Name' in descriptors sould match 'Molecule ChEMBL ID' in source
    # We allow for some format differences (e.g. trimming)
    
    mismatches = 0
    for idx in range(min(len(df_desc), len(df_source))):
        desc_name = str(df_desc.iloc[idx]['Name']).strip()
        source_id = str(df_source.iloc[idx]['Molecule ChEMBL ID']).strip()
        
        if desc_name != source_id:
            if mismatches < 5:
                print(f"Mismatch at index {idx}: Desc '{desc_name}' != Source '{source_id}'")
            mismatches += 1
            
    if mismatches == 0:
        print("SUCCESS: All IDs align perfectly.")
    else:
        print(f"FAILURE: {mismatches} mismatches found.")

    # Check if we can reproduce X and Y generation
    # The notebook drops 'Name' from X and takes 'pIC50' from a loaded df_padel
    # We need to see if we can calculate pIC50 from source or if we need to load it
    
    # Source has 'Standard Value' but usually pIC50 needs calculation
    # Let's inspect 'Standard Value'
    print("\nSample Standard Values:")
    print(df_source['Standard Value'].head())
    
    # Check if pIC50 is easily calculable (usually -log10(IC50_molar))
    # Standard Value is usually in nM.
    # pIC50 = -log10(Value * 10^-9) = 9 - log10(Value)
    
    # Let's try to match with the values seen in the notebook or the other csv
    if os.path.exists('Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv'):
        df_fp_pic50 = pd.read_csv('Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv')
        print(f"\nFinal dataset shape: {df_fp_pic50.shape}")
        if 'pIC50' in df_fp_pic50.columns:
            print("pIC50 column found in final dataset.")
            print(df_fp_pic50['pIC50'].head())
            
            # Check correlation if counts match
            if len(df_fp_pic50) == len(df_source):
                # We can't map exactly without IDs in df_fp_pic50, but we can assume order
                pass
        else:
            print("pIC50 column NOT found in final dataset?")

if __name__ == "__main__":
    check_alignment()
