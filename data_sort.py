import pandas as pd
import openpyxl

# Loading in the SNP data and the variant types that we consider non-coding
SNP_file = pd.read_excel('Data/Liu. et al.2023 studys 320 SNPs.xlsx',skiprows=1)
with open('Data/non_coding_types.txt','r') as f:
    non_coding_types= f.read().splitlines()

# Splitting into non-coding SNPs and coding SNPs
def is_non_coding(consequence):
    if pd.isna(consequence):
        return False
    consequence_list = [item.strip() for item in str(consequence).split(',')]
    return any(item in non_coding_types for item in consequence_list)

non_coding_mask = SNP_file['Consequence'].apply(is_non_coding)

non_coding_SNPs = SNP_file[non_coding_mask]
coding_SNPs = SNP_file[~non_coding_mask]

All_SNPS = []
CD_SNPS = []
UC_SNPS = []

SNP_file_CD = non_coding_SNPs[non_coding_SNPs['Phenotype_loci']=='CD']
SNP_file_UC = non_coding_SNPs[non_coding_SNPs['Phenotype_loci']=='UC']
SNP_file_IBD = non_coding_SNPs[non_coding_SNPs['Phenotype_loci']=='IBD']

for row in SNP_file_CD.itertuples():
    SNP_string = f"{row.RS_ID}_chr{row.Chr_index}_{row.Pos_index}_{row.A1_index}_{row.A2_index}"
    All_SNPS.append(SNP_string)
    CD_SNPS.append(SNP_string)

for row in SNP_file_UC.itertuples():
    SNP_string = f"{row.RS_ID}_chr{row.Chr_index}_{row.Pos_index}_{row.A1_index}_{row.A2_index}"
    All_SNPS.append(SNP_string)
    UC_SNPS.append(SNP_string)

for row in SNP_file_IBD.itertuples():
    SNP_string = f"{row.RS_ID}_chr{row.Chr_index}_{row.Pos_index}_{row.A1_index}_{row.A2_index}"
    All_SNPS.append(SNP_string)
    UC_SNPS.append(SNP_string)
    CD_SNPS.append(SNP_string)

All_SNPS_unique = list(dict.fromkeys(All_SNPS))
CD_SNPS_unique = list(dict.fromkeys(CD_SNPS))
UC_SNPS_unique = list(dict.fromkeys(UC_SNPS))

with open('Data/IBD_SNPS.txt','w') as file:
    for item in All_SNPS_unique:
        file.write(item + "\n")

with open('Data/CD_SNPS.txt','w') as file:
    for item in CD_SNPS_unique:
        file.write(item + "\n")

with open('Data/UC_SNPS.txt','w') as file:
    for item in UC_SNPS_unique:
        file.write(item + "\n")




