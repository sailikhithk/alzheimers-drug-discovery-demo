# Linux/Shell Commands Used in Project

## File Operations

| Command | Purpose | Usage in Notebook |
|---------|---------|-------------------|
| `cat` | Display file contents | `! cat molecule.smi \| head -5` — Preview SMILES file |
| `head` | Show first N lines | `! cat molecule.smi \| head -5` — Show first 5 molecules |
| `wc -l` | Count lines in file | `! cat molecule.smi \| wc -l` — Count molecules |
| `ls -l` | List files with details | `! ls -l` — Verify descriptors_output.csv exists |
| `chmod +x` | Make file executable | `! chmod +x Miniconda3-py37_4.8.2-Linux-x86_64.sh` |
| `unzip` | Extract zip archive | `! unzip padel.zip` — Extract PaDEL descriptor tool |

## Package Management

| Command | Purpose | Usage in Notebook |
|---------|---------|-------------------|
| `pip install` | Install Python packages | `! pip install chembl_webresource_client` |
| `pip uninstall -y` | Remove packages without prompt | `! pip uninstall -y numpy scikit-learn` |
| `conda install` | Install via Conda | `! conda install -c rdkit rdkit -y` |

## Download & Network

| Command | Purpose | Usage in Notebook |
|---------|---------|-------------------|
| `wget` | Download files from URL | `! wget https://github.com/.../padel.zip` |

## Script Execution

| Command | Purpose | Usage in Notebook |
|---------|---------|-------------------|
| `bash` | Run shell script | `! bash ./Miniconda3-py37_4.8.2-Linux-x86_64.sh -b -f -p /usr/local` |
| `bash padel.sh` | Run PaDEL fingerprint generation | Generates `descriptors_output.csv` from `molecule.smi` |

## Full Commands Reference

### Package Installation (Commented in main notebook)
```bash
! pip install "numpy<2.0.0" scikit-learn==1.2.2
! pip install lazypredict
! pip install chembl_webresource_client
! pip install rdkit
! pip install openpyxl
! pip install padel-pywrapper
! pip install padelpy
```

### Miniconda/RDKit Installation (Kaggle-specific, commented out)
```bash
! wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh
! chmod +x Miniconda3-py37_4.8.2-Linux-x86_64.sh
! bash ./Miniconda3-py37_4.8.2-Linux-x86_64.sh -b -f -p /usr/local
! conda install -c rdkit rdkit -y
```

### PaDEL Descriptor Setup (Kaggle-specific, commented out)
```bash
! wget https://github.com/dataprofessor/bioinformatics/raw/master/padel.zip
! wget https://github.com/dataprofessor/bioinformatics/raw/master/padel.sh
! unzip padel.zip
```

### Data Verification Commands
```bash
! cat molecule.smi | head -5      # Preview first 5 SMILES
! cat molecule.smi | wc -l        # Count total molecules
! ls -l                           # Verify output files exist
```

## Notes

- Commands prefixed with `!` run in shell from Jupyter notebook cells
- Many commands are commented out in the main notebook (for local execution)
- Kaggle/Colab environments may require different installation approaches
- PaDEL requires Java runtime to be installed
- For local macOS execution, use `pip install rdkit` instead of conda
