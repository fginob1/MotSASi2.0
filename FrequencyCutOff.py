import pandas as pd
from Bio import SeqIO
import statistics
import json
import math
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
from itertools import repeat

# obtenemos el nombre del gen
def gene_names(xml_record, sp_h_genes_primary_name_path):
    """
    This function takes the UniProt SeqRecord object and returns the codifying gene primary name as also the
    documented synonyms, as long as they does not match with the gene primary name of another human protein
    in SwissProt
    """

    with open(sp_h_genes_primary_name_path, "r") as file:
        primary_name_genes = json.load(file)

    genes = []

    if "gene_name_primary" in xml_record.annotations.keys():
        genes.append(xml_record.annotations["gene_name_primary"])
    if "gene_name_synonym" in xml_record.annotations.keys():
        for gene in xml_record.annotations["gene_name_synonym"]:
            if gene not in primary_name_genes:
                genes.append(gene)

    return genes

def cut_off_freq(UniProtID, xml_db_path, sp_h_genes_primary_name_path, cv_her_af_path):
    """
    This function takes the UniProtID of the query protein and returns the -log of the allelic frequency
    of the pathogenic variant in ClinVar with such evidence that makes it difficult to clearly stablish
    its pathogenicity. Said this, it is a threshold from which we cannot assure that if the allelic frequency
    is greater, its ClinVar reports will probably point the variant as pathogenic for mendelian disorders
    If there are too few variants, we return a frequency cut-off bi heritage model or a general ClinVar cut-off
    """
    print(UniProtID)

    # load the UniProt, SwissProt, human, database
    xmls = SeqIO.index_db(f"{xml_db_path}.idx", f"{xml_db_path}", "uniprot-xml")

    # load the dataframe with ClinVar variants annotated with frequency and heritage models
    cv_her_af = pd.read_csv(cv_her_af_path, dtype={"Cromosoma": str})
    #print(cv_her_af)

    # load the UniProt SeqRecord and get the gene names of the protein
    xml_record = xmls[UniProtID]
    gene = gene_names(xml_record, sp_h_genes_primary_name_path)

    # filter variants by gene names
    gene_variants = cv_her_af[cv_her_af.Gene.isin(gene)]

    switch = False

    # if there are no variants, return a general ClinVar frequency cut-off
    if gene_variants.empty:
        freq_cutoff = 4.1
        switch = True

    # if there are more than 10 pathogenic variants, we will iterate over them
    elif len(gene_variants[gene_variants.Simple_ClinSig == "Pathogenic"]) >= 10:

        # filter by pathogenic variants
        gene_variants = gene_variants[gene_variants.Simple_ClinSig == "Pathogenic"]
        
        # sort variants by its descending frequency
        gene_variants = gene_variants.sort_values("AF_popmax", ascending=False)

        # iterate through filtered variants
        for i, variant in gene_variants.iterrows():

            # very curated variants are clearly pathogenic, continue to the next variant
            if variant.ClinSig in ["reviewed by expert panel", "practice guideline", "association"]:
                continue
            else:

                # if the variant has more than 25% of benign references, it gives us a lot of doubts about its
                # pathogenicity, continue to the next variant
                if variant.Benign_references/variant.Total_references >= 0.25:
                    continue
                else:

                    # if the variant has more than 5 references, we can continue with the next variant
                    if variant.Total_references >= 5:
                        continue
                    else:

                        # if the variant has more than 50% of uncertain references, it gives us a lot of doubts
                        # about its pathogenicity, continue to the next variant
                        if variant.Uncertain_references/variant.Total_references >= 0.5:
                            continue
                        else:

                            # we ask the variant to have 50% or more of pathogenic references, if not continue
                            if variant.Pathogenic_references/variant.Total_references < 0.5:
                                continue
                            else:
                                
                                # check if the allele frequency is a float number (not a np.nan cell value)
                                if not np.isnan(variant.AF_popmax):
                                    # if the variant pass the previous filters, return the -log of the allelic frequency
                                    freq_cutoff = -math.log10(variant.AF_popmax)
                                    switch = True
                                    break

    if switch == False:

        # there are variants but few, lets ask the heritage model of the gene
        inheritance = gene_variants.iloc[0].inheritance

        # assign a freq cut-off for each heritage model, if the gene has not a stablished one, we use a
        # general ClinVar frequency cut-off
        if inheritance == "AR":
            freq_cutoff = 4.1
            pass
        elif inheritance == "AD":
            freq_cutoff = 4.28
            pass
        elif inheritance == "AD/AR":
            freq_cutoff = 4.18
            pass
        elif np.isnan(inheritance):
            freq_cutoff = 4.1

    return [UniProtID, freq_cutoff]

def generate_genes_dictionary(xml_db_path, proteoma_list_path, output_json_path):

    xmls = SeqIO.index_db(f"{xml_db_path}.idx", f"{xml_db_path}", "uniprot-xml")

    raw_sp = pd.read_csv(proteoma_list_path, header=None)
    lista_sp = list(set(raw_sp[0]))
    print(f"Cantidad de UniProt IDs únicos: {len(lista_sp)}")

    genes = []

    for UniProtID in tqdm(lista_sp):
        xml_record = xmls[UniProtID]
        if "gene_name_primary" in xml_record.annotations:
            genes.append(xml_record.annotations["gene_name_primary"])

    with open(output_json_path, "w") as file:
        json.dump(genes, file)

    return genes


if __name__ == "__main__":

    xml_db_path = "./UniProt/uniprot_sprot_h.xml"
    sp_h_genes_primary_name_path = "./UniProt/uniprot_primary_name_genes.json"
    cv_her_af_path = "./ClinVar/clinvar_mix.csv"    
    proteoma_list_path = "./UniProt/uniprot_sprot_h.list"

    generate_genes_dictionary(xml_db_path, proteoma_list_path, sp_h_genes_primary_name_path)

    proteoma_sp = pd.read_csv(f"{proteoma_list_path}", header=None)
    lista_proteoma_sp = list(proteoma_sp[0])

    thr_dict = {}

    with Pool() as pool:
        cut_off = pool.starmap(cut_off_freq, zip(lista_proteoma_sp, repeat(xml_db_path), repeat(sp_h_genes_primary_name_path), repeat(cv_her_af_path)))

    df = pd.DataFrame(cut_off, columns=["UniProtID", "Freq_cut_off"])
    df.to_csv("./GnomAD/uniprot_genes_freq_cutoffs.csv", index=False)
