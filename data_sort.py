import pandas as pd
import openpyxl

SNP_file_1 = pd.read_excel('Data/167 IBD SNPs.xlsx',sheet_name = 'Detailed assoc stats', skiprows= range(164,176))
SNP_file_2 = pd.read_excel('Data/241 GWAS IBD SNPs.xlsx',skiprows=8)
SNP_file_3 = pd.read_excel('Data/Liu. et al.2023 studys 320 SNPs.xlsx',skiprows=1)
SNP_file_4 = pd.read_excel('Data/scIBD database 318 risk SNPs.xlsx')

All_SNPS = []
CD_SNPS = []
UC_SNPS = []

# #  SNP file 1, To do: add chromosome position, cycle through UC and IBD and add to respective lists
# SNP_file_1_CD = SNP_file_1[SNP_file_1['type']=='CD']
# SNP_file_1_UC = SNP_file_1[SNP_file_1['type']=='UC']
# SNP_file_1_IBD = SNP_file_1[SNP_file_1['type']=='IBD']

# for row in SNP_file_1_CD.itertuples():
#     SNP_string = f"{row.GWAS_SNP}_chr{row.Chr}_{row.GWAS_nonrisk}_{row.GWAS_risk}"
#     All_SNPS.append(SNP_string)
#     UC_SNPS.append(SNP_string)

# SNP file 2
    
# SNP file 3

SNP_file_3_CD = SNP_file_3[SNP_file_3['Phenotype_loci']=='CD']
SNP_file_3_UC = SNP_file_3[SNP_file_3['Phenotype_loci']=='UC']
SNP_file_3_IBD = SNP_file_3[SNP_file_3['Phenotype_loci']=='IBD']

for row in SNP_file_3_CD.itertuples():
    SNP_string = f"{row.RS_ID}_chr{row.Chr_index}_{row.Pos_index}_{row.A1_index}_{row.A2_index}"
    All_SNPS.append(SNP_string)
    CD_SNPS.append(SNP_string)

for row in SNP_file_3_UC.itertuples():
    SNP_string = f"{row.RS_ID}_chr{row.Chr_index}_{row.Pos_index}_{row.A1_index}_{row.A2_index}"
    All_SNPS.append(SNP_string)
    UC_SNPS.append(SNP_string)

for row in SNP_file_3_IBD.itertuples():
    SNP_string = f"{row.RS_ID}_chr{row.Chr_index}_{row.Pos_index}_{row.A1_index}_{row.A2_index}"
    All_SNPS.append(SNP_string)
    UC_SNPS.append(SNP_string)
    CD_SNPS.append(SNP_string)

# Writing to text files and remocing duplicates

All_SNPS_unique = list(dict.fromkeys(All_SNPS))
CD_SNPS_unique = list(dict.fromkeys(CD_SNPS))
UC_SNPS_unique = list(dict.fromkeys(UC_SNPS))


with open('Data/IBD_SNPS.txt','w') as file:
    for item in All_SNPS:
        file.write(item + "\n")

with open('Data/CD_SNPS.txt','w') as file:
    for item in CD_SNPS:
        file.write(item + "\n")

with open('Data/UC_SNPS.txt','w') as file:
    for item in UC_SNPS:
        file.write(item + "\n")




