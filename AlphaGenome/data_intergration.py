"""
This program takes the outputs of SNP_list_data_load.py and website_data_load.py in the format RSID_CHR_POS_REF_ALT_TYPE.
This program merges the two files together so that the new file contains saved as merged_SNP_dataset:
- A list of RSIDs from the raw SNP file with the type data added from the files from the website. This is only matched based
 on the RSID and not neccessarily the nucleotide change
- SNPs with no disease type from both files are left blank
- SNPs with multiple different phenotypes associated have all these phenotypes in the last column

To do:
- Validate
- Prevent duplicates i.e. CD,CD


"""
import pandas as pd
import numpy as np

# Loads the data from each file
website_SNPs = pd.read_csv('AlphaGenome/Results/dataset_combination/website_SNP_list.txt',
                            sep='_', header=None, 
                            names=['RSID','CHR', 'POS','REF','ALT','TYPE'])
Raw_SNPs = pd.read_csv('AlphaGenome/Results/dataset_combination/raw_SNP_list.txt',
                        sep='_', header=None, 
                        names=['RSID','CHR','POS','REF','ALT','TYPE'])

# Reduces the website SNP dataframe to just RSID and type as these are the only relevant data to add to Raw_SNPs file
website_SNP_refined = website_SNPs[['RSID', 'TYPE']]
website_SNP_refined= website_SNP_refined.drop_duplicates()
website_SNP_refined.to_csv('AlphaGenome/Results/test_todelete/intergration_test.txt',index=None,header=None)
# Combines the type column in the refined SNP dataset i.e CD,UC
# website_SNP_refined= website_SNP_refined.groupby('RSID')['TYPE'].apply(
#     lambda x: ','.join(x.dropna().unique())
# ).reset_index()

# website_SNP_refined.to_csv('AlphaGenome/Results/test_todelete/intergration_test.txt',index=None,header=None)
# website_SNP_refined['TYPE'] = website_SNP_refined['TYPE'].replace('', np.nan)

# Raw_SNPs= Raw_SNPs.drop_duplicates()

# merged_SNPs = pd.merge(Raw_SNPs, website_SNP_refined, on='RSID', how='left', suffixes=('_raw','_web'))

# def combine_types(row):
#     raw_val = row['TYPE_raw']
#     web_val = row['TYPE_web']
    
#     if pd.notna(raw_val) and pd.notna(web_val):
#         return str(raw_val) + ',' + str(web_val)
#     elif pd.isna(raw_val) and pd.notna(web_val):
#         return web_val
#     elif pd.notna(raw_val) and pd.isna(web_val):
#         return raw_val
#     else:  # both are NaN
#         return np.nan

# merged_SNPs['Type'] = merged_SNPs.apply(combine_types, axis=1)

# merged_SNPs = merged_SNPs.drop(columns=['TYPE_raw', 'TYPE_web'])
# merged_SNPs = merged_SNPs.drop_duplicates()

# merged_SNP_string = merged_SNPs.astype(str).agg('_'.join, axis=1)

# print(f"Count of NaN in type: {merged_SNPs['Type'].isna().sum()}")
# with open('AlphaGenome/Results/dataset_combination/merged_SNP_dataset.txt', 'w') as f:
#     for line in merged_SNP_string:
#         f.write(line + '\n')


