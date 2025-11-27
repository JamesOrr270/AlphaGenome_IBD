"""
This program creates lists of patient specific SNPs based on a matrix of SNPs with the functions:
- Takes only the patients from the matrix that we have gene expression data for
- For each patient converts their matrix of SNPs into a text file of their SNPs to be used in SNP_affected_genes.py

Issues:
- Only am getting 39 patients should be getting 48
- Found that in the summary spreadsheet some are repeated so there are only 41 unique VICCD
- Also found that two in the summary spread spreadsheet do not appear in the SNP dataset [9875, 7908]
"""

import pandas as pd


patient_genetic_data = pd.read_csv('Gene_expression/genetic_data/SNPs_uc_gandalf_revision.txt', sep='\t')
patient_summary = pd.read_csv('Gene_expression/expression_data/Microarray_Leuven_UC_uninflamed.txt', sep='\t')


# Getting the list of patients we have genetic data for

patient_list = []

for row in patient_summary.itertuples():
    patient_list.append(int(row.VLECCID))

# Getting the SNPs for each patient list and writing them to a text file


for patient, SNP in patient_genetic_data.items():

    if '_' not in patient:
        continue

    if int(patient.split('_')[1]) in patient_list:
        patient_SNPs = patient_genetic_data.loc[SNP[SNP == 1].index,'SNP'].tolist()

        with open(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_SNP_lists/{patient.split('_')[1]}_SNP_list.txt','w') as f:
            for snp in patient_SNPs:
                f.write(f"{snp}\n")


