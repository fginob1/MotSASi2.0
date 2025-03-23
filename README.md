# MotSASi 2.0

MotSASi 2.0 is a computational pipeline designed to enhance the prediction of functional Short Linear Motifs (SLiMs) by integrating genomic variant information and structural data from both PDB-deposited crystallographic structures and AlphaFold2 models. Additionally, it allows the generation of Single Amino Acid Substitution (SAS) matrices, providing predictions for all missense variants occurring within SLiMs.

If you find this tool useful, we kindly ask that you cite our work.

---

## Requirements and Installation

For better compatibility and reproducibility, we recommend running the scripts within an Anaconda environment.

### Install Anaconda
```bash
wget https://repo.anaconda.com/archive/Anaconda3-2024.02-0-Linux-x86_64.sh
bash Anaconda3-2024.02-0-Linux-x86_64.sh
```

### Create and Activate the Environment
```bash
conda create --name motsasi python=3.10
conda activate motsasi
```

### Install Required Dependencies
```bash
conda install -c conda-forge biopython
```

### Clone the Repository
```bash
git clone https://github.com/fginob1/MotSASi2.0.git
cd MotSASi2.0
git lfs pull
```

### Create Necessary Directories
```bash
mkdir tmp Motifs PDB AF_human_proteome repaired_pdbs secondary_structures scratch
```

---

## Required Software

### FoldX
To run FoldX locally, request an academic license at: [FoldX Academic License](https://foldxsuite.crg.eu/academic-license-info)

Ensure the **FoldX executable** and the **"molecules" folder** are in the main `MotSASi2.0` directory.

### FreeSASA
```bash
pip install freesasa
```

### SCRATCH
```bash
cd scratch
wget http://download.igb.uci.edu/SCRATCH-1D_1.3.tar.gz
tar -zxf SCRATCH-1D_1.3.tar.gz
cd SCRATCH-1D_1.3
perl install.pl
cd ../..
```

---

## Required Datasets

### PDB Database (Zipped)
**CAUTION:** If you already have the zipped PDB database (`pdbxxxx.ent.gz` files) on your computer, specify the correct path instead of redownloading.

```bash
cd PDB
mkdir zipped unzipped
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/all/pdb/ ./zipped
cd ..
```

### Human Proteome AlphaFold Database
```bash
cd AF_human_proteome
wget -c https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
tar -xvf UP000005640_9606_HUMAN_v4.tar
rm UP000005640_9606_HUMAN_v4.tar
cd ..
```

### SwissProt and Trembl FASTA Files from UniProt
```bash
cd UniProt
wget -O uniprot_sprot.fasta.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
zcat uniprot_sprot.fasta.gz | bgzip > uniprot_sprot.fasta.bgz

wget -O uniprot_trembl.fasta.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
zcat uniprot_trembl.fasta.gz | bgzip > uniprot_trembl.fasta.bgz
cd ..
```

---

## Running MotSASi 2.0

### `MotSASi_2_cluster.py`: Integrating Structural & Variant Information for SLiM Analysis

The script `MotSASi_2_cluster.py` performs several key tasks:

1. **Motif Identification**: Searches for the motif of interest within the human proteome.
2. **Variant Collection**: Retrieves benign and pathogenic variants from ClinVar and gnomAD.
3. **Substitution Matrix Calculation**:
   - If PDB structures are available, the matrix is generated using FoldX.
   - If not, an AlphaFold2 model prediction is used.
4. **ClinSig and Allele Frequency Matrix Generation**: Based on ELM motif instances (used as a positive control).
5. **Conservation Score Calculation** for the positive control.
6. **Secondary Structure Probability Calculation** using SCRATCH for the positive control.
7. **Solvent Accessible Surface Area (SASA) Calculation** for the positive control.
8. **GO Terms Collection** from UniProt for the positive control.

### Running the Script
To execute `MotSASi_2_cluster.py`, run the following command:

```bash
python3 MotSASi_2_cluster.py [MOTIF] [DOT-SEPARATED MOTIF] [ELM MOTIF NAME]
```

#### Example Command
```bash
python3 MotSASi_2_cluster.py [RK]P[^P][^P]L.[LIVMF] RK.P.^P.^P.L.x.LIVMF DOC_MAPK_JIP1_4
```

Where:
- `[MOTIF]` → The motif in regular expression format.
- `[DOT-SEPARATED MOTIF]` → The motif with elements separated by dots.
- `[ELM MOTIF NAME]` → The name of the motif in the ELM database.

### Output
Once executed, the script generates a motif-specific folder containing all output files. This folder is stored in the parent directory of the folder where the script is run.

---

This README provides a structured and methodical guide for setting up and running the MotSASi 2.0 pipeline. If you have any issues, please reach out via GitHub issues.

