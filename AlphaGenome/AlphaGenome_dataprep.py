import pandas as pd

with open('AlphaGenome/Data/SNP_rawdata_yufan/GWAS_study_GRCh38.txt','r') as file:
    SNP1 = file.read().strip().split('\n')

with open('AlphaGenome/Data/SNP_rawdata_yufan/IBD_SNPs_list_GRCh38.txt','r') as file:
    SNP2 = file.read().strip().split('\n')

with open('AlphaGenome/Data/SNP_rawdata_yufan/Liu_et_al_2023_GRCh38.txt','r') as file:
    SNP3 = file.read().strip().split('\n')

with open('AlphaGenome/Data/SNP_rawdata_yufan/Original_iSNP_CD_list_GRCh38.txt','r') as file:
    SNP4 = file.read().strip().split('\n')

with open('AlphaGenome/Data/SNP_rawdata_yufan/Original_iSNP_UC_list_GRCh38.txt','r') as file:
    SNP5 = file.read().strip().split('\n')

combined_SNPs = SNP1 + SNP2 + SNP3 + SNP4 + SNP5
total_SNPs = pd.Series(combined_SNPs).drop_duplicates()

with open('AlphaGenome/Results/AlphaGenome_dataprep/totalSNPs.txt','w') as file:
    for SNP in total_SNPs:
        file.write(f"{SNP}\n")
