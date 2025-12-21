
import json
import os

notebook_path = 'ml-for-alzheimer-s-drug-discovery-in-progress.ipynb'

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

changed = False

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        new_source = []
        for line in cell['source']:
            original_line = line
            # Specific replacements based on grep results
            if "'/kaggle/input/beta-amyloid-a4-protein-active-compounds/" in line:
                 line = line.replace('/kaggle/input/beta-amyloid-a4-protein-active-compounds/', '')
            
            if "'/kaggle/input'" in line:
                line = line.replace("'/kaggle/input'", "'.'")
            
            if line != original_line:
                changed = True
                print(f"Modifying line: {original_line.strip()} -> {line.strip()}")
            new_source.append(line)
        cell['source'] = new_source

if changed:
    print("Writing changes to notebook...")
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(nb, f, indent=1)
    print("Done.")
else:
    print("No changes needed.")
