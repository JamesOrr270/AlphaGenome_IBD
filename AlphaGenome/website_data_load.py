"""
This program reads the data from the four excel files and writes them to a txt file as a list of strings in the format
RSID_CHR_POS_REF_ALT_DIS. If the information is not presant in the excel file it is intended to return NaN for that column

It functions through the create_snp_series function which takes the arguments:
df - the datagrame of the excel file
RSID - Name of the column in the excel file containing RSID
CHR - Name of the column in the excel file containing the chromosome number
POS - Name of the column in the excel file containing the position of the SNP
REF - Name of the column containing the reference variant
ALT - Name of th ecolumn containing the alternate variant
TYPE - Name of the column containing the diseses type i.e. IBD, CD or UC
"""

import pandas as pd

def create_snp_df(df, RSID=None, CHR=None, POS=None, REF=None, ALT=None, TYPE=None):

    SNP_df = pd.DataFrame(columns = ['RSID','CHR','POS','REF','ALT','TYPE'])

    if RSID is not None:
        SNP_df['RSID'] = df[RSID]
    
    if CHR is not None:
        SNP_df['CHR'] = "CHR" + df[CHR].astype(str)

    if POS is not None:
        SNP_df['POS'] = df[POS]
    
    if REF is not None:
        SNP_df['REF'] = df[REF]
    
    if ALT is not None:
        SNP_df['ALT'] = df[ALT]
    
    if TYPE is not None:
        SNP_df['TYPE'] = df[TYPE]

    SNP_df['TYPE'] = SNP_df['TYPE'].str.replace(',',' ')
    SNP_df['RSID'] = SNP_df['RSID'].str.replace('†', '').str.replace('§', '')
    
    return SNP_df

#########################################################################################################################
# Liu et al 2023
SNP_file_1 = pd.read_excel('AlphaGenome/Data/SNP_Raw_Data_from_Website/Liu. et al.2023 studys 320 SNPs.xlsx',skiprows=1)

SNP_list_1 = create_snp_df(SNP_file_1,RSID='RS_ID',CHR='Chr_index',POS='Pos_index',REF='A1_index',ALT='A2_index',TYPE='Phenotype_loci')

#########################################################################################################################
# 241 IBD SNPS

SNP_file_2 = pd.read_excel('AlphaGenome/Data/SNP_Raw_Data_from_Website/241 GWAS IBD SNPs.xlsx',skiprows=8)

SNP_list_2 = create_snp_df(SNP_file_2,RSID='topSNP rsid',CHR='Chr',POS='topSNP Position (bp)',TYPE='Trait')

#########################################################################################################################
# 167 IBD SNPS

SNP_file_3 = pd.read_excel('AlphaGenome/Data/SNP_Raw_Data_from_Website/167 IBD SNPs.xlsx',sheet_name= 'Detailed assoc stats',skipfooter=11)

SNP_list_3 = create_snp_df(SNP_file_3,RSID='GWAS_SNP',CHR='Chr',REF='GWAS_nonrisk',ALT='GWAS_risk',TYPE='type')

#########################################################################################################################
# Asian SNPs

SNP_file_4 = pd.read_excel('AlphaGenome/Data/SNP_Raw_Data_from_Website/Asian SNPs.xlsx')

SNP_list_4 = create_snp_df(SNP_file_4,RSID='SNP',CHR='Chr',POS='Position',TYPE='Type')

#########################################################################################################################
# Combined list

combined_SNPs= pd.concat([SNP_list_1,SNP_list_2,SNP_list_3,SNP_list_4])
combined_SNPs.to_csv('AlphaGenome/Results/dataset_combination/website_SNP_list.txt',sep='_',index=False,header=False)
