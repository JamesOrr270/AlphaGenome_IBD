"""
This program takes the text file containing the SNPs from the raw files and anoather text file containing the SNPs from 
the website both in the format RSID_CHR_POS_REF_ALT_TYPE. Merges the two files together so that the new file contains:
- A list of RSIDs from the raw SNP file with the type data added from the files from the website. This is only matched based
 on the RSID and not neccessarily the nucleotide change
- SNPs with no disease type from both files are left blank
- SNPs with multiple different phenotypes associated have all these phenotypes in the last column

To do:
- Validate
- Prevent duplicates i.e. CD,CD
- Find out why there are more SNPs in the merged dataset compared to the Raw SNPs list

"""
import pandas as pd
import numpy as np

website_SNPs = pd.read_csv('AlphaGenome/Results/dataset_combination/website_SNP_list.txt',
                            sep='_', header=None, 
                            names=['RSID','CHR', 'POS','REF','ALT','TYPE'])
Raw_SNPs = pd.read_csv('AlphaGenome/Results/dataset_combination/SNP_list_Yufan.txt',
                        sep='_', header=None, 
                        names=['RSID','CHR','POS','REF','ALT','TYPE'])

website_disease = website_SNPs[['RSID', 'TYPE']]
website_disease = website_disease.drop_duplicates()


website_disease = website_SNPs.groupby('RSID')['TYPE'].apply(
    lambda x: ','.join(x.dropna().unique())
).reset_index()

website_disease['TYPE'] = website_disease['TYPE'].replace('', np.nan)

Raw_SNPs= Raw_SNPs.drop_duplicates()

merged_SNPs = pd.merge(Raw_SNPs, website_disease, on='RSID', how='left', suffixes=('_raw','_web'))

def combine_types(row):
    raw_val = row['TYPE_raw']
    web_val = row['TYPE_web']
    
    if pd.notna(raw_val) and pd.notna(web_val):
        return str(raw_val) + ',' + str(web_val)
    elif pd.isna(raw_val) and pd.notna(web_val):
        return web_val
    elif pd.notna(raw_val) and pd.isna(web_val):
        return raw_val
    else:  # both are NaN
        return np.nan

merged_SNPs['Type'] = merged_SNPs.apply(combine_types, axis=1)

merged_SNPs = merged_SNPs.drop(columns=['TYPE_raw', 'TYPE_web'])
merged_SNPs = merged_SNPs.drop_duplicates()

merged_SNP_string = merged_SNPs.astype(str).agg('_'.join, axis=1)

# print(f"Count of NaN in type: {merged_SNPs['Type'].isna().sum()}")
with open('AlphaGenome/Results/dataset_combination/merged_SNP_dataset.txt', 'w') as f:
    for line in merged_SNP_string:
        f.write(line + '\n')


