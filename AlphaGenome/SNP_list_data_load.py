"""
This file takes the text files that yufan sent me and and outs them all in the same file writing them to a new text file.
It also adds an extra term at the end which is blank for most the files but for the UC and CD files adds this on the end

Issues:
-For some reason the amount of SNPs has increased 
"""

import pandas as pd

with open('AlphaGenome/Data/SNP_rawdata_yufan/GWAS_study_GRCh38.txt','r') as file:
    SNP1 = [line + "_" for line in file.read().strip().split('\n')]

with open('AlphaGenome/Data/SNP_rawdata_yufan/IBD_SNPs_list_GRCh38.txt','r') as file:
    SNP2 = [line + "_" for line in file.read().strip().split('\n')]

with open('AlphaGenome/Data/SNP_rawdata_yufan/Liu_et_al_2023_GRCh38.txt','r') as file:
    SNP3 = [line + "_" for line in file.read().strip().split('\n')]

with open('AlphaGenome/Data/SNP_rawdata_yufan/Original_iSNP_CD_list_GRCh38.txt','r') as file:
    SNP4 = [line + "_CD" for line in file.read().strip().split('\n')]

with open('AlphaGenome/Data/SNP_rawdata_yufan/Original_iSNP_UC_list_GRCh38.txt','r') as file:
    SNP5 = [line + "_UC" for line in file.read().strip().split('\n')]

combined_SNPs = SNP1 + SNP2 + SNP3 + SNP4 + SNP5
print(f"Combined SNPs: {len(combined_SNPs)}")

total_SNPs = pd.Series(combined_SNPs).drop_duplicates()

print(f"Total SNPs: {len(total_SNPs)}")
with open('AlphaGenome/Results/dataset_combination/SNP_list_Yufan.txt','w') as file:
    for SNP in total_SNPs:
        file.write(f"{SNP}\n")

