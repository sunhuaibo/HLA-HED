#!/usr/bin/env python
# -*- coding=utf-8 -*-

# =====================================
# Author: Huaibo Sun
# E-mail: huaibo_sun@foxmail.com
# date:   2022-03-31
# =====================================

import pandas as pd
from Bio import SeqIO
from pathlib import Path
from collections import Counter
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def get_opt():
    """
    Input HLA file format
    
    Sample A1      A2      B1      B2      C1      C2
    p1     A*01:01 A*01:03 B*07:01 B*07:02 C*01:01 C*01:02
    p2     A*01:01 A*01:03 B*07:01 B*07:02 C*01:01 C*01:02
    
    If you use this tool, please cite the following three papers.
    
    Grantham R. Amino acid difference formula to help explain protein evolution. Science. 1974 Sep 6;185(4154):862-4. doi: 10.1126/science.185.4154.862. PMID: 4843792.
    Pierini F, Lenz TL. Divergent Allele Advantage at Human MHC Genes: Signatures of Past and Ongoing Selection. Mol Biol Evol. 2018 Sep 1;35(9):2145-2158. doi: 10.1093/molbev/msy116. PMID: 29893875; PMCID: PMC6106954.
    Chowell D, Krishna C, Pierini F, Makarov V, Rizvi NA, Kuo F, Morris LGT, Riaz N, Lenz TL, Chan TA. Evolutionary divergence of HLA class I genotype impacts efficacy of cancer immunotherapy. Nat Med. 2019 Nov;25(11):1715-1720. doi: 10.1038/s41591-019-0639-4. Epub 2019 Nov 7. PMID: 31700181; PMCID: PMC7938381.
    
    """
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, epilog=get_opt.__doc__)
    parser.add_argument("-d", required=True, help="Distance matrix for all amino acids (reference: DOI: 10.1126/science.185.4154.862)")
    parser.add_argument("-f", required=True, help="Amino acid sequences in fasta format")
    parser.add_argument("-i", required=True, help="Input file of tab-delimited with individual HLA typing")
    parser.add_argument("-o", required=True, help="Output file name")

    parse = parser.parse_args()
    return(parse)

def check_file(infile):
    if not infile.exists:
        raise(f"{str(infile)} file is not exist")

def read_fasta(infile):
    infile = Path(infile)
    check_file(infile)
    record = SeqIO.parse(infile, "fasta")
    seq_array = {seq.id: str(seq.seq) for seq in record}
    seq_len = [len(value) for value in seq_array.values()]
    if len(Counter(seq_len)) != 1:
        raise("Input sequences length is not equality")
    return(seq_array)

def read_aa(infile):
    infile = Path(infile)
    check_file(infile)
    df = pd.read_csv(infile, header=0, sep="\t", index_col=0)
    aa_pairwise_dis = df.to_dict()
    return(aa_pairwise_dis)

def calculate_distange(hla1, hla2, sequence, distance):
    seq_hla1 = sequence.get(hla1, False)
    seq_hla2 = sequence.get(hla2, False)
    seq_len = len(seq_hla1)
    if not seq_hla1 or not seq_hla2:
        return("NA")
    else:
        dis = 0
        for i in range(seq_len):
            aa1 = seq_hla1[i]
            aa2 = seq_hla2[i]
            dis += distance[aa1][aa2]
        dis = dis / seq_len
        return(dis)


def main():
    opt = get_opt()
    seq_array = read_fasta(opt.f)
    aa_pairwise_dis = read_aa(opt.d)

    infile = Path(opt.i)
    outfile = Path(opt.o)
    check_file(infile)

    df = pd.read_csv(infile, header=0, sep="\t")
    
    outheader = ["Sample","HED_A","HED_B","HED_C","Mean_HE"]
    with open(outfile, "w") as fw:
        fw.write("\t".join(outheader) + "\n")
        for _, line in df.iterrows():
            hla_a1 = line["A1"]
            hla_a2 = line["A2"]
            dis_hla_a = calculate_distange(hla_a1, hla_a2, seq_array, aa_pairwise_dis)

            hla_b1 = line["B1"]
            hla_b2 = line["B2"]
            dis_hla_b = calculate_distange(hla_b1, hla_b2, seq_array, aa_pairwise_dis)
            
            hla_c1 = line["C1"]
            hla_c2 = line["C2"]
            dis_hla_c = calculate_distange(hla_c1, hla_c2, seq_array, aa_pairwise_dis)

            if dis_hla_a == "NA" or dis_hla_b == "NA" or dis_hla_c == "NA":
                dis_mean = "NA"
            else:
                dis_mean = (dis_hla_a + dis_hla_b + dis_hla_c) / 3

            outline = [line["Sample"], dis_hla_a, dis_hla_b, dis_hla_c, dis_mean]
            outline = [str(x) for x in outline]

            fw.write("\t".join(outline) + "\n")

if __name__ == "__main__":
    main()