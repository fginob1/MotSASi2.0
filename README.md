# MotSASi 2.0

This repository contains all the necessary code files to run the MotSASi 2.0 pipeline.

MotSASi 2.0 is a strategy designed to enhance the prediction of functional Short Linear Motifs (SLiMs) by integrating genomic variant information and structural data from both PDB-deposited crystallographic structures and AlphaFold2 models. Additionally, it enables the generation of Single Amino Acid Substitution (SAS) matrices, providing predictions for all missense variants occurring within SLiMs.

If you find this tool useful, we kindly ask you to cite our work.

---

## Requirements and Setup

For better compatibility and reproducibility, we recommend running the scripts within an Anaconda environment.

### **Install Anaconda and Create Environment**

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2024.02-0-Linux-x86_64.sh
bash Anaconda3-2024.02-0-Linux-x86_64.sh
conda create --name motsasi python=3.10
conda activate motsasi
```

### **Install Dependencies**

Biopython is required:

```bash
conda install -c conda-forge biopython mafft
```

### **Clone the Repository**

```bash
git clone https://github.com/fginob1/MotSASi2.0.git
cd MotSASi2.0
git lfs pull
```

### **Create Required Directories**

```bash
mkdir tmp PDB AF_human_proteome repaired_pdbs secondary_structures scratch
```

### **Install Required Software**

#### **FoldX**

Go to [FoldX Website](https://foldxsuite.crg.eu/academic-license-info) and request an academic license to run FoldX locally. Important files that must be present in the main directory `MotSASi2.0` are the FoldX executable and the "molecules" folder.

#### **FreeSASA**

```bash
pip install freesasa
```

#### **ssbio**

```bash
pip install ssbio
```


#### **SCRATCH**

```bash
cd scratch
wget http://download.igb.uci.edu/SCRATCH-1D_1.3.tar.gz
tar -zxf SCRATCH-1D_1.3.tar.gz
cd SCRATCH-1D_1.3
perl install.pl
cd ../..
```

There is an issue with Biopython in order to run ssbio to parse scratch output. This is represented by an impossibility to import IUPAC from Bio.Alphabet. To resolve this issue, we simply should remove the imports and usages of IUPAC in the utils.py file in the ssbio folder. We provide three Linux terminal lines in order to change this file with no functional alterations. The user must specify the specific path to the folder (if working with an environment, the folder should probably be there).

```bash
sed -i '/from Bio.Alphabet import IUPAC/d' ~/path_to_the_ssbbio_folder/ssbio/protein/sequence/utils/utils.py

sed -i 's/def cast_to_seq(obj, alphabet=IUPAC.extended_protein):/def cast_to_seq(obj):/' ~/path_to_the_ssbbio_folder/ssbio/protein/sequence/utils/utils.py

sed -i 's/def cast_to_seq_record(obj, alphabet=IUPAC.extended_protein,/def cast_to_seq_record(obj,/' ~/path_to_the_ssbbio_folder/ssbio/protein/sequence/utils/utils.py
```

---

## Required Datasets

### **Zipped PDB Database**

**CAUTION**: If you already have the zipped PDB database (`pdbxxxx.ent.gz` files) on your computer, just specify the correct path. The same applies to the Human proteome AlphaFold database and the SwissProt/Trembl files from UniProt.

```bash
mkdir -p PDB/zipped PDB/unzipped PDB/divided
rsync -avP rsync.wwpdb.org::ftp/data/structures/divided/pdb/ "./PDB/divided/"
find "./PDB/divided/" -type f -name "*.ent.gz" -exec mv -n {} "./PDB/zipped/" \;
rm -r ./PDB/divided
```

### **Human Proteome AlphaFold Database**

```bash
cd AF_human_proteome
wget -c https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
tar -xvf UP000005640_9606_HUMAN_v4.tar
rm UP000005640_9606_HUMAN_v4.tar
cd ..
```

### **SwissProt and Trembl FASTA Files from UniProt**

```bash
cd UniProt
wget -O "uniprot_sprot.fasta.gz" https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
zcat uniprot_sprot.fasta.gz | bgzip > uniprot_sprot.fasta.bgz
wget -O "uniprot_trembl.fasta.gz" https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
zcat uniprot_trembl.fasta.gz | bgzip > uniprot_trembl.fasta.bgz
wget -O uniprot_sprot.list "https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28*%29+AND+%28reviewed%3Atrue%29"
wget -O uniprot_human_sprot.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28*%29+AND+%28model_organism%3A9606%29+AND+%28reviewed%3Atrue%29"
wget -O uniprot_human_sprot.list "https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28*%29+AND+%28model_organism%3A9606%29+AND+%28reviewed%3Atrue%29"
wget -O uniprot_human_sprot.xml "https://rest.uniprot.org/uniprotkb/stream?format=xml&query=%28*%29+AND+%28model_organism%3A9606%29+AND+%28reviewed%3Atrue%29"
cd ..
```

### **Decompress Required Files**

To extract the necessary `.gz` files, run the following command:

```bash
gunzip ./ClinVar/ClinVar_missense_all_filtered.csv.gz \
       ./GnomAD/GnomAD_missense.csv.gz \
       ./Predictions/motsasi_variant_predictions.csv.gz
```

---

## Running MotSASi 2.0

### **MotSASi_2_cluster.py: Integrating Structural & Variant Information for SLiM Analysis**

This script performs the following tasks:

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

### **Execution**

```bash
python3 MotSASi_2.py [ELM MOTIF NAME] [MOTIF] [DOT-SEPARATED MOTIF]
```

**Example:**

```bash
python3 MotSASi_2.py DOC_MAPK_JIP1_4 [RK]P[^P][^P]L.[LIVMF] RK.P.^P.^P.L.x.LIVMF
```

---

## Incorporating AlphaFold2 Models

When no crystallographic structure is available for a motif class, AlphaFold-generated models can be used. Users can generate these models via the [AlphaFold Server](https://alphafoldserver.com/) or install AlphaFold locally.

### **Adding AlphaFold Models**

Place `.pdb` prediction files inside a `Seed` folder within the corresponding motif directory. An example is given for the CLV_C14_Caspase3-7 motif class:
   ```bash
   ./Motifs/CLV_C14_Caspase3-7/Seed
   ```
### **Running AlphaFold Parameter Extraction**

Execute the following command:

```bash
python3 alphafold_parameters.py [ELM MOTIF NAME]
```

**Example:**

```bash
python3 alphafold_parameters.py CLV_C14_Caspase3-7
```

This generates a TSV file named 'interface_res_stat_total.tsv' containing Interaction Energy and Confidence metrics used to select optimal models. Once completed, you can run the pipeline with AlphaFold2 models integrated into the FoldX matrices.
