"""
This program reads the data from the four excel files and writes them to a txt file as a list of strings in the format
RSID_CHR_POS_REF_ALT_DIS. 

It functions through the create_snp_series function which takes the arguments:
df - the datagrame of the excel file
RSID - Name of the column in the excel file containing RSID
CHR - Name of the column in the excel file containing the chromosome number
POS - Name of the column in the excel file containing the position of the SNP
REF - Name of the column containing the reference variant
ALT - Name of th ecolumn containing the alternate variant
disease - Name of the column containing the diseses type i.e. IBD, CD or UC
"""

import pandas as pd

def create_snp_series(df, RSID, CHR, POS, REF, ALT, disease):
    # Create a list to store formatted SNP strings
    snp_list = []
    
    for idx, row in df.iterrows():
        # Handle 'None' string by converting to empty string
        rsid = ' ' if RSID == 'None' else str(row[RSID])
        chr_val = ' ' if CHR == 'None' else str(row[CHR])
        pos = ' ' if POS == 'None' else str(row[POS])
        ref = ' ' if REF == 'None' else str(row[REF])
        alt = ' ' if ALT == 'None' else str(row[ALT])
        dis = ' ' if disease == 'None' else str(row[disease])
        
        snp_string = f'{rsid}_chr{chr_val}_{pos}_{ref}_{alt}_{dis}'
        snp_list.append(snp_string)
    
    return pd.Series(snp_list)

#########################################################################################################################
# Liu et al 2023
SNP_file_1 = pd.read_excel('AlphaGenome/Data/SNP_Raw_Data_from_Website/Liu. et al.2023 studys 320 SNPs.xlsx',skiprows=1)

SNP_list_1 = create_snp_series(SNP_file_1,'RS_ID','Chr_index','Pos_index','A1_index','A2_index','Phenotype_loci')

SNP_list_1.to_csv('AlphaGenome/Results/test_todelete/SNP_list_1.txt',index=False, header=False)


#########################################################################################################################
# 241 IBD SNPS

SNP_file_2 = pd.read_excel('AlphaGenome/Data/SNP_Raw_Data_from_Website/241 GWAS IBD SNPs.xlsx',skiprows=8)

SNP_list_2 = create_snp_series(SNP_file_2,'topSNP rsid','Chr','topSNP Position (bp)','None','None','Trait')

SNP_list_2.to_csv('AlphaGenome/Results/test_todelete/SNP_list_2.txt',index=False, header=False)

#########################################################################################################################
# 167 IBD SNPS

SNP_file_3 = pd.read_excel('AlphaGenome/Data/SNP_Raw_Data_from_Website/167 IBD SNPs.xlsx',sheet_name= 'Detailed assoc stats',skipfooter=11)

SNP_list_3 = create_snp_series(SNP_file_3,'GWAS_SNP','Chr','None','GWAS_nonrisk','GWAS_risk','type')

SNP_list_3.to_csv('AlphaGenome/Results/test_todelete/SNP_list_3.txt',index=False, header=False)

#########################################################################################################################
# Asian SNPs

SNP_file_4 = pd.read_excel('AlphaGenome/Data/SNP_Raw_Data_from_Website/Asian SNPs.xlsx')

SNP_list_4 = create_snp_series(SNP_file_4,'SNP','Chr','Position','None','None','Type')

SNP_list_4.to_csv('AlphaGenome/Results/test_todelete/SNP_list_4.txt',index=False, header=False)

#########################################################################################################################
# Combined list

combined_SNP_list = pd.concat([SNP_list_1,SNP_list_2,SNP_list_3,SNP_list_4])
combined_SNP_list.to_csv('AlphaGenome/Results/dataset_combination/website_SNP_list.txt',index=False,header=False)
