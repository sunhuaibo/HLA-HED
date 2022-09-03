#!/usr/bin/env python
# -*- coding=utf-8 -*-

# =====================================
# Author: Huaibo Sun
# E-mail: huaibo_sun@foxmail.com
# date:   2022-03-31
# =====================================

import os
from typing import NoReturn
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from itertools import combinations
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter

def get_opt() -> Namespace:
    """
    Input HLA file format
    
    Sample A1      A2      B1      B2      C1      C2
    p1     A*01:01 A*01:03 B*07:01 B*07:02 C*01:01 C*01:02
    p2     A*01:01 A*01:03 B*07:01 B*07:02 C*01:01 C*01:02
    
    If you use this tool, please cite the three following papers.
    
    Grantham R. Amino acid difference formula to help explain protein evolution. Science. 1974 Sep 6;185(4154):862-4. doi: 10.1126/science.185.4154.862. PMID: 4843792.
    Pierini F, Lenz TL. Divergent Allele Advantage at Human MHC Genes: Signatures of Past and Ongoing Selection. Mol Biol Evol. 2018 Sep 1;35(9):2145-2158. doi: 10.1093/molbev/msy116. PMID: 29893875; PMCID: PMC6106954.
    Chowell D, Krishna C, Pierini F, Makarov V, Rizvi NA, Kuo F, Morris LGT, Riaz N, Lenz TL, Chan TA. Evolutionary divergence of HLA class I genotype impacts efficacy of cancer immunotherapy. Nat Med. 2019 Nov;25(11):1715-1720. doi: 10.1038/s41591-019-0639-4. Epub 2019 Nov 7. PMID: 31700181; PMCID: PMC7938381.
    
    """
    
    script = os.path.dirname(os.path.abspath(__file__))
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, epilog=get_opt.__doc__)
    parser.add_argument("-d", default=f"{script}/database/grantham_matrix.txt", help="Distance matrix for all amino acids, default: database/grantham_matrix.txt. (reference: DOI: 10.1126/science.185.4154.862)")
    parser.add_argument("-f", default=f"{script}/database/ABC_prot.fa", help="Amino acid sequences in fasta format, default: database/ABC_prot.fa.")
    parser.add_argument("-i", required=True, help="Input file of tab-delimited with individual HLA typing.")
    parser.add_argument("-p", action="store_true", help="Paired HED score.")
    parser.add_argument("-o", required=True, help="Output file name.")

    parse = parser.parse_args()
    return parse

def check_file(infile: str) -> NoReturn:
    if not infile.exists:
        raise Exception(f"{str(infile)} file does not exist")

def read_fasta(infile: str) -> dict:
    infile = Path(infile)
    check_file(infile)
    record = SeqIO.parse(infile, "fasta")
    sequences = {seq.id: str(seq.seq) for seq in record}
    seq_len = [len(value) for value in sequences.values()]
    if len(set(seq_len)) != 1:
        raise Exception("Input sequences length is not equality")
    return sequences

def read_aa(infile: str) -> dict:
    infile = Path(infile)
    check_file(infile)
    df = pd.read_csv(infile, header=0, sep="\t", index_col=0)
    aa_pairwise_dis = df.to_dict()
    return aa_pairwise_dis

def calculate_distange(hla1: str, hla2: str, sequences: dict, distance: dict) -> float:
    seq_hla1 = sequences.get(hla1, False)
    seq_hla2 = sequences.get(hla2, False)
    if not seq_hla1 or not seq_hla2:
        return("NA")
    else:
        seq_len = len(seq_hla1)
        dis = 0
        for i in range(seq_len):
            aa1 = seq_hla1[i]
            aa2 = seq_hla2[i]
            dis += distance[aa1][aa2]
        dis = dis / seq_len
        return dis


def main():
    opt = get_opt()
    sequences = read_fasta(opt.f)
    aa_pairwise_dis = read_aa(opt.d)

    infile = Path(opt.i)
    outfile = Path(opt.o)
    check_file(infile)

    df = pd.read_csv(infile, header=0, sep="\t")

    if opt.p:
        df2 = pd.melt(df, id_vars=["Sample"], value_vars=["A1", "A2", "B1","B2", "C1","C2"])
        alleles = set(df2["value"].values.tolist())
        alleles_pair = combinations(alleles, 2)
    
        outheader = ["Allele1","Allele2","HED"]
        with open(outfile, "w") as fw:
            fw.write("\t".join(outheader) + "\n")
            for allele1, allele2 in alleles_pair:
                dis_hla_pair = calculate_distange(allele1, allele2, sequences, aa_pairwise_dis)
                outline = [allele1, allele2, dis_hla_pair]
                outline = [str(x) for x in outline]

                fw.write("\t".join(outline) + "\n")
    else:
        outheader = ["Sample","HED_A","HED_B","HED_C","Mean_HE"]
        with open(outfile, "w") as fw:
            fw.write("\t".join(outheader) + "\n")
            for _, line in df.iterrows():
                hla_a1 = line["A1"]
                hla_a2 = line["A2"]
                dis_hla_a = calculate_distange(hla_a1, hla_a2, sequences, aa_pairwise_dis)

                hla_b1 = line["B1"]
                hla_b2 = line["B2"]
                dis_hla_b = calculate_distange(hla_b1, hla_b2, sequences, aa_pairwise_dis)
                
                hla_c1 = line["C1"]
                hla_c2 = line["C2"]
                dis_hla_c = calculate_distange(hla_c1, hla_c2, sequences, aa_pairwise_dis)

                if dis_hla_a == "NA" or dis_hla_b == "NA" or dis_hla_c == "NA":
                    dis_mean = "NA"
                else:
                    dis_mean = (dis_hla_a + dis_hla_b + dis_hla_c) / 3

                outline = [line["Sample"], dis_hla_a, dis_hla_b, dis_hla_c, dis_mean]
                outline = [str(x) for x in outline]

                fw.write("\t".join(outline) + "\n")

if __name__ == "__main__":
    main()