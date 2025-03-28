from Bio.PDB import PDBParser
import pandas as pd
import statistics
import argparse
import re
import os

def annotate_parameters(motif_name, combinations_df, repaired_pdbs_path, cwd_path, tmp_path):
    """
    This script takes the previously generated csv file and starts iterating over the AlfhaFold models
    and calculating the required parameters, interaction energy and confidence of the motif residues
    """

    combinations = []

    for n, row in combinations_df.iterrows():

        # INTERACTION ENERGY

        # check if the crystal structure has already been repaired, prevent losing time
        if not os.path.exists(f"{repaired_pdbs_path}{row.file_name}_Repair.pdb"):
            # FoldX can only make calculations if the crystal structure file is in the working directory
            if not os.path.exists(f"{cwd_path}{row.file_name}.pdb"):
                os.system(f"cp {cwd_path}Motifs/{motif_name}/Seed/{row.file_name}.pdb {cwd_path}")
            print(f"Repairing {row.file} AlphaFold model")
            # repair the crystal structure, outputs all the files in tmp and then we move it to the working directory
            os.system(f"./foldx --command=RepairPDB --pdb={row.file_name}.pdb --output-dir={tmp_path} > {tmp_path}output_foldx.txt")
            os.system(f"mv {tmp_path}{row.file_name}_Repair.pdb {cwd_path}")
            os.system(f"rm {tmp_path}{row.file_name}_Repair.fxout")

        else:
            os.system(f"mv {repaired_pdbs_path}{row.file_name}_Repair.pdb {cwd_path}")
            os.system(f"cp {cwd_path}Motifs/{motif_name}/Seed/{row.file_name}.pdb {cwd_path}")

        # calculate the binding ddG
        os.system(f"./foldx --command=AnalyseComplex --pdb={row.file_name}_Repair.pdb --analyseComplexChains=A,B --output-dir={tmp_path} > {tmp_path}output_foldx.txt")
        # we dont need anymore the repaired pdb file, we move it with the others repaired pdbs
        if not os.path.exists(f"{repaired_pdbs_path}{row.file_name}_Repair.pdb"):
            os.system(f"mv {row.file_name}_Repair.pdb {repaired_pdbs_path}")
        ac_df = pd.read_csv(f"{tmp_path}Summary_{row.file_name}_Repair_AC.fxout", sep="\t", skiprows=8)
        interaction_energy = ac_df["Interaction Energy"][0]

        # RESIDUES CONFIDENCE

        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure("alphafold", f"{cwd_path}{row.file_name}.pdb")

        motif_residues_confidence = []

        for model in structure:
            for chain in model:
                if chain.id == "B":
                    for residue in chain:
                        for atom in residue:
                            if atom.name == "CA":
                                motif_residues_confidence.append(atom.bfactor)

        mean_motif_residue_confidence = statistics.mean(motif_residues_confidence)

        os.system(f"rm {cwd_path}{row.file_name}.pdb")

        combinations.append((motif_name, row.rank, row.model, interaction_energy, mean_motif_residue_confidence))

    combinations_df = pd.DataFrame(combinations, columns=["Motif_name", "rank", "model", "Interaction_energy", "DDT_mean_motif_IR"])
    combinations_df.to_csv(f"{cwd_path}Motifs/{motif_name}/interface_res_stat_total.csv")

    return combinations_df

def initialize_csv_file(motif_name, cwd_path):
    """
    This scripts creates a csv file from the AlphaFold models in the Seed folder
    """

    combinations = []

    for file in os.listdir(f"{cwd_path}Motifs/{motif_name}/Seed/"):
        if file.endswith(".pdb"):
            match = re.search(r"rank_(\d+)_model_(\d+)", file)
            rank_number = int(match.group(1))
            model_number = int(match.group(2))
            combinations.append((motif_name, rank_number, model_number, file[:-4]))

    combinations_df = pd.DataFrame(combinations, columns=["Motif_name", "rank", "model", "file_name"])
    combinations_df.to_csv(f"{cwd_path}Motifs/{motif_name}/interface_res_stat_raw.csv")

    return combinations_df

def alphafold_parameters(motif_name, repaired_pdbs_path, cwd_path, tmp_path):
    """
    This scripts firt creats a raw csv only with rankings and models, then it fills it with the required metrics
    in order to filter AlphaFold models
    """

    combinations_df = initialize_csv_file(motif_name, cwd_path)
    combinations_df = annotate_parameters(motif_name, combinations_df, repaired_pdbs_path, cwd_path, tmp_path)

    return combinations_df

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Annotate a CSV file with some model parameters in order to filter bad AlphaFold models and retain good AlphaFold models")
    parser.add_argument("Motif_name", action="store", help = "Name of the motif, e.g. CLV_C14_Caspase3-7")
    parser.add_argument("-rep_pdbs", "--repaired_pdbs", action="store", help = "folder with repaired PDBSs", default="./repaired_pdbs/")
    parser.add_argument("-cwd", "--cwd", action="store", help = "current working directory path", default="./")
    parser.add_argument("-tmp", "--tmp_path", action="store", help = "folder receiving foldx output", default="./tmp/")

    args = parser.parse_args()

    _ = alphafold_parameters(args.Motif_name, args.repaired_pdbs, args.cwd, args.tmp_path)
