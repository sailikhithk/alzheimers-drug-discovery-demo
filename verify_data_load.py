
import pandas as pd
import os

file_path = 'Beta_amyloid A4_protein_active_compounds.xlsx'

print(f"Checking if {file_path} exists...")
if os.path.exists(file_path):
    print("File exists.")
    try:
        df = pd.read_excel(file_path)
        print(f"Successfully loaded {file_path}")
        print(f"Data shape: {df.shape}")
        print("First 5 rows:")
        print(df.head())
    except Exception as e:
        print(f"Error loading file: {e}")
else:
    print("File not found.")

print("\nChecking os.walk('.') output:")
for dirname, _, filenames in os.walk('.'):
    for filename in filenames:
        if filename.endswith('.csv') or filename.endswith('.xlsx'):
            print(os.path.join(dirname, filename))
