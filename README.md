# MotSASi 2.0

This repository contains all the necessary code files to run the MotSASi 2.0 pipeline.

MotSASi 2.0 is a strategy designed to enhance the prediction of functional Short Linear Motifs (SLiMs) by integrating genomic variant information and structural data from both PDB-deposited crystallographic structures and AlphaFold2 models. Additionally, it enables the generation of Single Amino Acid Substitution (SAS) matrices, providing predictions for all missense variants occurring within SLiMs.

If you find this tool useful, we kindly ask you to cite our work.

---

## Requirements and Usage Instructions

For better compatibility and reproducibility, we recommend running the scripts within an Anaconda environment.

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2024.02-0-Linux-x86_64.sh
bash Anaconda3-2024.02-0-Linux-x86_64.sh
conda create --name motsasi python=3.10
conda activate motsasi
```

### Required Dependencies

Biopython is required:

```bash
conda install -c conda-forge biopython
```

### Clone the Repository

```bash
git clone https://github.com/fginob1/MotSASi2.0.git
cd MotSASi2.0
git lfs pull
```

### Create Required Directories

```bash
mkdir tmp Motifs PDB AF_human_proteome repaired_pdbs secondary_structures scratch
```

### Install Required Software

#### FoldX

Go to [FoldX Website](https://foldxsuite.crg.eu/academic-license-info) and request an academic license to run FoldX locally. Important files that must be present in the main directory `MotSASi2.0` are the FoldX executable and the "molecules" folder.

#### FreeSASA

```bash
pip install freesasa
```

#### SCRATCH

```bash
cd scratch
wget http://download.igb.uci.edu/SCRATCH-1D_1.3.tar.gz
tar -zxf SCRATCH-1D_1.3.tar.gz
cd SCRATCH-1D_1.3
perl install.pl
cd ../..
```

### Required Datasets

#### Zipped PDB Database

**CAUTION**: If you already have the zipped PDB database (`pdbxxxx.ent.gz` files) on your computer, just specify the correct path. The same applies to the Human proteome AlphaFold database and the SwissProt/Trembl files from UniProt.

```bash
cd PDB
mkdir zipped unzipped
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/all/pdb/ ./zipped
cd ..
```

#### Human Proteome AlphaFold Database

```bash
cd AF_human_proteome
wget -c https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
tar -xvf UP000005640_9606_HUMAN_v4.tar
rm UP000005640_9606_HUMAN_v4.tar
cd ..
```

#### SwissProt and Trembl FASTA Files from UniProt

```bash
cd UniProt
wget -O "uniprot_sprot.fasta.gz" https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
zcat uniprot_sprot.fasta.gz | bgzip > uniprot_sprot.fasta.bgz
wget -O "uniprot_trembl.fasta.gz" https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
zcat uniprot_trembl.fasta.gz | bgzip > uniprot_trembl.fasta.bgz
cd ..
```

---

## Running MotSASi 2.0

### MotSASi_2_cluster.py: Integrating Structural & Variant Information for SLiM Analysis

The script `MotSASi_2_cluster.py` performs the following tasks:

- **Motif Identification**: Searches for the motif of interest within the human proteome.
- **Variant Collection**: Retrieves benign and pathogenic variants from ClinVar and gnomAD.
- **Substitution Matrix Calculation**:
  - If PDB structures are available, the matrix is generated using FoldX.
  - If not, an AlphaFold2 model prediction is used.
- **ClinSig Matrix and Allele Frequency Matrix Generation**: Based on ELM motif instances (used as a positive control).
- **Conservation Score Calculation** for the positive control.
- **Secondary Structure Probability Calculation** using SCRATCH for the positive control.
- **Solvent Accessible Surface Area (SASA) Calculation** for the positive control.
- **GO Terms Collection** from UniProt for the positive control.

### Execution

Run the following command:

```bash
python3 MotSASi_2.py [MOTIF] [DOT-SEPARATED MOTIF] [ELM MOTIF NAME]
```

**Example Command:**

```bash
python3 MotSASi_2.py "[RK]P[^P][^P]L.[LIVMF]" "RK.P.^P.^P.L.x.LIVMF" "DOC_MAPK_JIP1_4"
```

Where:
- `[MOTIF]` → The motif in regular expression format.
- `[DOT-SEPARATED MOTIF]` → The motif with elements separated by dots.
- `[ELM MOTIF NAME]` → The name of the motif in the ELM database.

### Output

Once executed, the script generates a motif-specific folder containing all output files. This folder is stored in the parent directory of the folder where the script is run.

---

## Introducing AlphaFold2 Models into the Pipeline

When a crystallographic structure is not associated with a given motif class, AlphaFold-generated models can be incorporated into the pipeline. The user can generate these models via the [AlphaFold Server](https://alphafoldserver.com/) or by installing AlphaFold locally. According to the MotSASi 2.0 methodology, users should run AlphaFold using the sequence of a given motif along with the sequence of its binding domain.

### Adding AlphaFold Models

Once the user has the `.pdb` files with predictions, they should be placed in a folder named `Seed` inside the corresponding motif folder within `Motifs`.

For example, for the **DOC_MAPK_JIP1_4** motif class:

```bash
./Motifs/DOC_MAPK_JIP1_4/Seed
```

Additionally, a TSV file must be prepared with the following column names:

```
Motif	rank	model
```

This TSV file should be placed inside the motif folder:

```bash
./Motifs/DOC_MAPK_JIP1_4/
```

### Running AlphaFold Parameter Extraction

From the main directory (`./`), execute:

```bash
python3 alphafold_parameters.py [ELM MOTIF NAME]
```

**Example Command:**

```bash
python3 alphafold_parameters.py DOC_MAPK_JIP1_4
```

This step is necessary to prepare the TSV file with Interaction Energy and Confidence metrics, which will be used to select the optimal models. An example file for the DOC_MAPK_JIP1_4 motif class is provided in the Motifs folder. After this step, you can run the complete pipeline, which will incorporate the AlphaFold2 models to build the FoldX matrices.
