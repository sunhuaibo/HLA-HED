
## HLA-HED
---

A python script to calculate the HLA Evolutionary Divergence (HED) Score.


> Briefly, The protein sequence of exons 2 and 3 of each allele of  HLA-I genotype was extracted; these sequences correspond to the peptide-binding domains. Divergences (HED score) between allele sequences were calculated using the Grantham distance metric.

### Python dependencies

If you haven't already satisfied these dependencies on your system, install these Python (python 3.x) packages via pip or conda:

- [pandas](https://pandas.pydata.org/)
- [Biopython](https://biopython.org/)

### Database

- ABC_prot.fa, download from [IMGT](https://www.ebi.ac.uk/ipd/imgt/hla/), update 2022-01-14
- grantham_matrix.txt, from [Amino acid difference formula to help explain protein evolution](https://www.science.org/doi/10.1126/science.185.4154.862)

### Uasge

```
python hla_hed.py -h

optional arguments:
  -h, --help  show this help message and exit
  -d D        Distance matrix for all amino acids (reference: DOI: 10.1126/science.185.4154.862)
  -f F        Amino acid sequences in fasta format
  -i I        Input file of tab-delimited with individual HLA typing
  -p          Paired HED score
  -o O        Output file name

    Input HLA file format

    Sample A1      A2      B1      B2      C1      C2
    p1     A*01:01 A*01:03 B*07:01 B*07:02 C*01:01 C*01:02
    p2     A*01:01 A*01:03 B*07:01 B*07:02 C*01:01 C*01:02

    If you use this tool, please cite the following three papers.

    Grantham R. Amino acid difference formula to help explain protein evolution. Science. 1974 Sep 6;185(4154):862-4. doi: 10.1126/science.185.4154.862. PMID: 4843792.
    Pierini F, Lenz TL. Divergent Allele Advantage at Human MHC Genes: Signatures of Past and Ongoing Selection. Mol Biol Evol. 2018 Sep 1;35(9):2145-2158. doi: 10.1093/molbev/msy116. PMID: 29893875; PMCID: PMC6106954.
    Chowell D, Krishna C, Pierini F, Makarov V, Rizvi NA, Kuo F, Morris LGT, Riaz N, Lenz TL, Chan TA. Evolutionary divergence of HLA class I genotype impacts efficacy of cancer immunotherapy. Nat Med. 2019 Nov;25(11):1715-1720. doi: 10.1038/s41591-019-0639-4. Epub 2019 Nov 7. PMID: 31700181; PMCID: PMC7938381.

```

### Example

```
python hla_hed.py -d database/grantham_matrix.txt -f database/ABC_prot.fa -i example/test.txt -o example/test_hed.txt
```

```
python hla_hed.py -d database/grantham_matrix.txt -f database/ABC_prot.fa -i example/test.txt -p -o example/test_hed_paired.txt
```

### Citing
Please cite this paper when using HLA-HED for your publications.
> Jiang T, Jin Q, Wang J, Wu F, Chen J, Chen G, Huang Y, Chen J, Cheng Y, Wang Q, Pan Y, Zhou J, Shi J, Xu X, Lin L, Zhang W, Zhang Y, Liu Y, Fang Y, Feng J, Wang Z, Hu S, Fang J, Shu Y, Cui J, Hu Y, Yao W, Li X, Lin X, Wang R, Wang Y, Shi W, Feng G, Ni J, Mao B, Ren D, Sun H, Zhang H, Chen L, Zhou C, Ren S. HLA-I evolutionary divergence confers response to PD-1 blockade plus chemotherapy in untreated advanced non-small-cell lung cancer. Clin Cancer Res. 2023 Jul 14:CCR-23-0604. doi: 10.1158/1078-0432.CCR-23-0604. Epub ahead of print. PMID: 37449971.

### References

1. Grantham R. Amino acid difference formula to help explain protein evolution. Science. 1974 Sep 6;185(4154):862-4. doi: 10.1126/science.185.4154.862. PMID: 4843792.
2. Pierini F, Lenz TL. Divergent Allele Advantage at Human MHC Genes: Signatures of Past and Ongoing Selection. Mol Biol Evol. 2018 Sep 1;35(9):2145-2158. doi: 10.1093/molbev/msy116. PMID: 29893875; PMCID: PMC6106954.
3. Chowell D, Krishna C, Pierini F, Makarov V, Rizvi NA, Kuo F, Morris LGT, Riaz N, Lenz TL, Chan TA. Evolutionary divergence of HLA class I genotype impacts efficacy of cancer immunotherapy. Nat Med. 2019 Nov;25(11):1715-1720. doi: 10.1038/s41591-019-0639-4. Epub 2019 Nov 7. PMID: 31700181; PMCID: PMC7938381.
