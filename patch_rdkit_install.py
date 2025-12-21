
import json

notebook_path = 'ml-for-alzheimer-s-drug-discovery-in-progress.ipynb'

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

changed = False

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source_str = "".join(cell['source'])
        if "Miniconda3-py37_4.8.2-Linux-x86_64.sh" in source_str:
            print("Found target cell. Replacing...")
            cell['source'] = [
                "# Installing RDKit via pip (Recommended for local macOS)\n",
                "!pip install rdkit\n",
                "\n",
                "# Note: The original cell attempted to install Linux Miniconda,\n",
                "# which is incompatible with macOS and risky in a local notebook.\n"
            ]
            changed = True

if changed:
    print("Writing changes to notebook...")
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(nb, f, indent=1)
    print("Done.")
else:
    print("Target cell not found or already patched.")
