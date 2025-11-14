"""
This file takes the text files that yufan sent me and and outs them all in the same file writing them to a new text file.
It also adds an extra term at the end which is blank for most the files but for the UC and CD files adds this on the end

Caveat: 
- Have to make sure that the text files with the disease information are always concatonated secound as the logic of this program
keeps the last term of a match. 
"""

import pandas as pd

# Load Files
SNP1 = pd.read_csv('AlphaGenome/Data/SNP_rawdata_yufan/GWAS_study_GRCh38.txt', sep='_', names = ['RSID','CHR','POS','REF','ALT','DIS'],header=None)
SNP2 = pd.read_csv('AlphaGenome/Data/SNP_rawdata_yufan/IBD_SNPs_list_GRCh38.txt', sep='_', names = ['RSID','CHR','POS','REF','ALT','DIS'],header=None)
SNP3 = pd.read_csv('AlphaGenome/Data/SNP_rawdata_yufan/Liu_et_al_2023_GRCh38.txt', sep='_', names = ['RSID','CHR','POS','REF','ALT','DIS'],header=None)
SNP4 = pd.read_csv('AlphaGenome/Data/SNP_rawdata_yufan/Original_iSNP_CD_list_GRCh38.txt', sep='_', names = ['RSID','CHR','POS','REF','ALT','DIS'],header=None)
SNP5 = pd.read_csv('AlphaGenome/Data/SNP_rawdata_yufan/Original_iSNP_UC_list_GRCh38.txt', sep='_', names = ['RSID','CHR','POS','REF','ALT','DIS'],header=None)

# Add disease annotation to the relavant files
SNP4['DIS'] = 'CD'
SNP5['DIS'] = 'UC'

# concatonate the files and drop duplicates based on everthing but disease
total_SNPs = pd.concat([SNP1,SNP2,SNP3,SNP4,SNP5]).drop_duplicates(subset= ['RSID','CHR','POS','REF','ALT'],keep='last')

# Write to a text file in the correct format
total_SNPs.to_csv('AlphaGenome/Results/dataset_combination/raw_SNP_list.txt',sep='_',header=None,index=None)

