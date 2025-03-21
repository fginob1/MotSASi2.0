# MotSASi2.0
This repository contains all the necessary code files to run the MotSASi2.0 pipeline.

This repository contains all the necessary code files to run the MotSASi 2.0 pipeline.

MotSASi 2.0 is a strategy designed to enhance the prediction of functional Short Linear Motifs (SLiMs) by integrating genomic variant information and structural data from both PDB-deposited crystallographic structures and AlphaFold2 models. Additionally, it allows the generation of Single Amino Acid Substitution (SAS) matrices, providing predictions for all missense variants occurring within SLiMs.

If you find this tool useful, we kindly ask that you cite our work.

Requirements and Usage Instructions
To run MotSASi 2.0, the following software and tools must be installed:

- Python 3
- Biopython package
- FoldX
- FreeSASA
- SCRATCH

Additionally, the following datasets must be downloaded:

- Zipped PDB database
- Human proteome AlphaFold database
- SwissProt and Trembl FASTA files from UniProt

For better compatibility and reproducibility, we recommend running the scripts within an Anaconda Environment.

MotSASi_2_cluster.py: Integrating Structural & Variant Information for SLiM Analysis

The script MotSASi_2_cluster.py performs several tasks:

- Motif Identification: Searches for the motif of interest within the human proteome.
- Variant Collection: Retrieves benign and pathogenic variants from ClinVar and gnomAD.
- Substitution Matrix Calculation:
  - If PDB structures are available, the matrix is generated using FoldX.
  - If not, an AlphaFold2 model prediction is used.
- ClinSig Matrix and Allele Frequency Matrix Generation: Based on ELM motif instances (used as a positive control).
- Conservation Score Calculation for the positive control.
- Secondary Structure Probability Calculation using SCRATCH for the positive control.
- Solvent Accessible Surface Area (SASA) Calculation for the positive control.
- GO Terms Collection from UniProt for the positive control.

Running MotSASi 2.0

To execute MotSASi_2_cluster.py, follow these steps:

1) Open a Linux terminal.

2) Navigate to the folder containing the script (referred to as "motsasi" in the code).

Run the following command:

python3 MotSASi_2_cluster.py [MOTIF] [DOT-SEPARATED MOTIF] [ELM MOTIF NAME]

Example Command

python3 MotSASi_2_cluster.py [RK]P[^P][^P]L.[LIVMF] RK.P.^P.^P.L.x.LIVMF DOC_MAPK_JIP1_4

Where:

[MOTIF] → The motif in regular expression format.
[DOT-SEPARATED MOTIF] → The motif with elements separated by dots.
[ELM MOTIF NAME] → The name of the motif in the ELM database.

Output

Once executed, the script generates a motif-specific folder containing all output files. This folder is stored in the parent directory of the folder where the script is run.
