"""
Pipeline for GSEA for total patient SNPS (copy into main)
    get_total_patient_snps()
    run_alphaGenome_for_total_patient_snps()
    run_GSEA_for_total_patient_SNPs

Pipeline for ORA for total patient SNPs
    get_total_patient_snps()
    run_alphaGenome_for_total_patient_snps()
    overrepresentation_anaysis_for_total_patient_SNP()
"""


import gseapy as gp
import pandas as pd
from pybiomart import Dataset
from patient_variant_effect_Alphagenome import get_alphaGenome_prediction


patient_list = ['366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']
file_sizes = ['16KB','100KB','500KB','1MB']

def prepare_GSEA_data(expression_data_path):
    
    df = pd.read_csv(expression_data_path)
    
    # Finds the average raw score of samples that match on everything onther than biosample name
    step1 = df.groupby(
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

    ranking_df = step3[['gene_name', 'raw_score']].copy()
    ranking_df = ranking_df.sort_values(by='raw_score', ascending=False)

    return(ranking_df)
    
def perform_GSEA(gene_set,df):

   
    GSEA_res = gp.prerank(rnk=df,
                          gene_sets=gene_set,
                          seed=1024,
                          min_size= 5,
                          max_size= 1000)

    out_list = []

    for term in GSEA_res.results:
        p = GSEA_res.results[term]['pval']
        fdr = GSEA_res.results[term]['fdr']
        nes = GSEA_res.results[term]['nes']
        es = GSEA_res.results[term]['es']
        gene = GSEA_res.results[term]['lead_genes']
        out_list.append([term, p, fdr, nes, es, gene])

    df_out = pd.DataFrame(out_list, columns = ['Term','p_value','fdr', 'nes', 'es','gene']).sort_values('fdr').reset_index(drop = True)

    return(df_out)

def get_total_patient_snps():
       ################################# Getting a list of all SNPs presant in the patients for input into AlphaGenome
    total_SNP_df = pd.DataFrame()
    patient_list = ['366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']

    for patient in patient_list:
        df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_SNP_AG_input/{patient}_AG_SNPs.txt',
                         sep='_', names= ['RSID','CHROM','POS','REF','ALT'],
                         header=0)
        total_SNP_df = pd.concat([total_SNP_df,df],ignore_index=True)

    total_SNP_df = total_SNP_df.drop_duplicates()

    with open('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_SNP_AG_input/All_patient_SNPs.txt','w') as f:
        for _, row in total_SNP_df.iterrows():
            f.write(f"{row['RSID']}_{row['CHROM']}_{row['POS']}_{row['REF']}_{row['ALT']}\n")

def run_alphaGenome_for_total_patient_snps():
    ################################### Getting AlphaGenome Predictions
     file_sizes = ['16KB','100KB','500KB','1MB']

     for size in file_sizes:
            get_alphaGenome_prediction('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_SNP_AG_input/All_patient_SNPs.txt',
                                        f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_alphagenome_results/All_patient_{size}',
                                        size)
            
def run_GSEA_for_total_patient_SNPs():
    file_sizes = ['16KB','100KB','500KB','1MB']

    for size in file_sizes:
         GSEA_data = prepare_GSEA_data(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_alphagenome_results/All_patient_{size}_non_significant.csv')
        #  print(GSEA_data)
        #  GSEA_data = prepare_GSEA_data(f'//Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_SNPS_{size}_all_scores.csv')
         GSEA_results = perform_GSEA('GO_Biological_Process_2025',GSEA_data)
         GSEA_results = GSEA_results.sort_values('fdr', ascending=True)
         significant_results = GSEA_results[GSEA_results['fdr'] < 0.25]

         GSEA_results.to_csv((f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/GSEA_results/GSEA_results_{size}.csv'))
         significant_results.to_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/GSEA_results/GSEA_results_{size}_significant.csv')

def overrepresentation_anaysis_for_total_patient_SNP():
    for size in file_sizes:
        df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_alphagenome_results/All_patient_{size}_non_significant.csv')
        df = df[df['quantile_score']>0.99]
        df = df['gene_name']
        df = df.drop_duplicates()
        glist = df.squeeze().str.strip().to_list()

        enr = gp.enrichr(gene_list=glist, # or "./tests/data/gene_list.txt",
                 gene_sets='GO_Biological_Process_2025',
                 organism='human', 
                 outdir=None, 
                )
        
        enr.results.to_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/ORA_results/{size}_patients')

        enr.results = enr.results[enr.results['Adjusted P-value']< 0.05]
        enr.results.to_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/ORA_results/{size}_patients_significant')


if __name__ == '__main__':

    # get_total_patient_snps()
    # run_alphaGenome_for_total_patient_snps()
    # run_GSEA_for_total_patient_SNPs()

    overrepresentation_anaysis_for_total_patient_SNP()


 