"""
This code takes the txt files prouced from patient SNP list and the CSV files produced from Alphagenome.py and produces output file for each patient which contains the SNPs of each patient and
their corresponding gene expression data based on the Alphagenome predictions

Problems:
- Different variations of the same SNP have slighlty different regulation so hard to know which one to use (Could take an average)
- Some genes are regulated by multiple SNPs and so need to find a way of combining these (add these together)
- Some of the genes have different expression in different tissues in the colon but the sample name is just colon (maybe average these)
- Some of the 16KB do not have any and some do not have many therefore would it be better to use 100kb instead or could be used in addition

To Do:
- Look into the most effective way to aggregate scores
- Validate the score aggregation

Interesting Observation:
- Based on the different sensitivities seems like quite different predictions
"""
import pandas as pd
from tqdm import tqdm


# Splitting the variant ID column so that CHR, POS, REF, ALT are seperate and making sure all in the same formate and comaprable
def get_gene_predictions(filename,dataset):
    AlphaGenome_results = pd.read_csv(filename)
    split_cols = AlphaGenome_results['variant_id'].str.split(":", expand=True)

    AlphaGenome_results['CHR'] = split_cols[0] 
    AlphaGenome_results['POS'] = split_cols[1].astype(int)

    alleles = split_cols[2].str.split(">", expand=True)
    AlphaGenome_results['REF'] = alleles[0]
    AlphaGenome_results['ALT'] = alleles[1]  

    dataset = dataset.replace('','-').fillna('-')
    AlphaGenome_results['REF'] = AlphaGenome_results['REF'].replace('','-').fillna('-')
    AlphaGenome_results['ALT'] = AlphaGenome_results['ALT'].replace('','-').fillna('-')

    # Labelling results with RSID


    results_with_rsid = pd.merge(
        AlphaGenome_results,
        dataset[['RSID','CHR','POS','REF','ALT','DIS']],
        on=['CHR', 'POS', 'REF', 'ALT'],
        how='left'
    )

    results_with_rsid.to_csv('Gene_expression/test/AlphaGenome_with_RSID.csv')
    return(results_with_rsid)

def create_patients_gene_expression(filename,gene_prediction):
    patient_SNP_list = pd.read_csv(filename,names=['RSID'])
    patient_SNP_list['RSID'] = patient_SNP_list['RSID'].str.lower()

    patient_gene_expression = pd.merge(
        gene_prediction,
        patient_SNP_list,
        on = 'RSID',
        how = 'inner'
    )

    return(patient_gene_expression)

def score_aggregation(patients_gene_expression):
    patients_gene_expression = patients_gene_expression[['RSID', 'POS', 'REF', 'ALT', 'gene_name', 'gene_id', 'gene_type', 'biosample_name', 'raw_score']]
    
    # Finds the average raw score of samples that match on everything onther than biosample name
    step1 = patients_gene_expression.groupby(
        ['RSID', 'POS', 'REF', 'ALT', 'gene_name', 'gene_id', 'gene_type']
    ).agg({
        'raw_score': 'mean', 
    }).reset_index()

    # Finds the average raw score of those that match everything other than the REF and ALT
    step2 = step1.groupby(
        ['RSID', 'POS', 'gene_name', 'gene_id', 'gene_type']
    ).agg({
        'raw_score': 'mean',  
    }).reset_index()

    # Sums all affects of SNPs onto one gene
    step3 = step2.groupby(
        ['gene_name']
    ).agg({
        'raw_score': 'sum',  
    }).reset_index()

    patient_gene_list = step3[['gene_name','raw_score']]
    return(patient_gene_list)


SNP_dataset = pd.read_csv('AlphaGenome/Results/dataset_combination/merged_SNP_dataset.txt',sep='_',names=['RSID','CHR','POS','REF','ALT','DIS'])

patient_list = ['001','366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']
file_sizes = ['16KB','100KB','500KB','1MB']

for patient in tqdm(patient_list):
    for size in file_sizes:
        alphagenome_prediction = get_gene_predictions(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPs_{size}_all_scores.csv',SNP_dataset)
        gene_expression_df = create_patients_gene_expression(f'Gene_expression/patient_SNP_lists/{patient}_SNP_list.txt',alphagenome_prediction)
        score_aggregation(gene_expression_df).to_csv(f'Gene_expression/patient_gene_data_AG/{patient}_{size}_gene_list.csv',index=None, header=True)



