"""
Aim: I want to check if there is overlap between the SNPs in the patient dataset and the SNPs that were significant as of 
AlphaGenome. First need to add RSID back to the results file

What it does now:
- Takes the matrix of SNPs and writes SNPs presant to a text file - Need to validate 

To do:
- Check the overlap of RSIDs between the patient SNPs and alphagenome SNPs

Issues:
- Only am getting 39 patients should be getting 48
- Found that in the summary spreadsheet some are repeated so there are only 41 unique VICCD
- Also found that two in the summary spread spreadsheet do not appear in the SNP dataset [9875, 7908]
"""

import pandas as pd

patient_genetic_data = pd.read_csv('Gene_expression/genetic_data/SNPs_uc_gandalf_revision.txt', sep='\t')
iSNP_results = pd.read_csv('Gene_expression/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision.txt',sep='\t')
AlphaGenome_results = pd.read_csv('AlphaGenome/Results/AlphaGenome/All_SNPS_1MB_16NOV25_1539_all_scores.csv')
SNP_dataset = pd.read_csv('AlphaGenome/Results/dataset_combination/merged_SNP_dataset.txt',sep='_',names=['RSID','CHR','POS','REF','ALT','DIS'])
patient_summary = pd.read_csv('Gene_expression/expression_data/Microarray_Leuven_UC_uninflamed.txt', sep='\t')

# Splitting the variant ID column so that CHR, POS, REF, ALT are seperate and making sure all in the same formate and comaprable

split_cols = AlphaGenome_results['variant_id'].str.split(":", expand=True)

AlphaGenome_results['CHR'] = split_cols[0] 
AlphaGenome_results['POS'] = split_cols[1].astype(int)  

alleles = split_cols[2].str.split(">", expand=True)
AlphaGenome_results['REF'] = alleles[0]
AlphaGenome_results['ALT'] = alleles[1]  

SNP_dataset = SNP_dataset.replace('','-').fillna('-')
AlphaGenome_results['REF'] = AlphaGenome_results['REF'].replace('','-').fillna('-')
AlphaGenome_results['ALT'] = AlphaGenome_results['ALT'].replace('','-').fillna('-')

# Labelling results with RSID
results_with_rsid = pd.merge(
    AlphaGenome_results,
    SNP_dataset[['RSID','CHR','POS','REF','ALT']],
    on=['CHR', 'POS', 'REF', 'ALT'],
    how='left'
)

# Getting the list of patients we have genetic data for

patient_list = []

for row in patient_summary.itertuples():
    patient_list.append(int(row.VLECCID))

# Getting the SNPs for each patient list and writing them to a text file

count = 0

patient_list_2 = []
for patient, SNP in patient_genetic_data.items():

    if '_' not in patient:
        continue

    if int(patient.split('_')[1]) in patient_list:
        patient_SNPs = patient_genetic_data.loc[SNP[SNP == 1].index,'SNP'].tolist()

        with open(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/test/{patient.split('_')[1]}_SNP_list.txt','w') as f:
            for snp in patient_SNPs:
                f.write(f"{snp}\n")


