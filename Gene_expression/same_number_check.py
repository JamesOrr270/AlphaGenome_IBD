import pandas as pd

patient_list = ['001','366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']
all_patient_genes = pd.DataFrame()


for patient in patient_list:
    patient_expression = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_gene_data_AG/{patient}_1MB_gene_list.csv')
    patient_expression['patient'] = patient
    all_patient_genes = pd.concat([all_patient_genes,patient_expression])

# Reset index after concatenation
all_patient_genes = all_patient_genes.reset_index(drop=True)

# Check if patients share the same raw_score for the same gene
gene_score_check = all_patient_genes.groupby('gene_name')['raw_score'].nunique()

# Genes where all patients have the same raw_score
genes_same_score = gene_score_check[gene_score_check == 1]
print(f"\nGenes where all patients share the same raw_score: {len(genes_same_score)}")

# Genes where patients have different raw_scores
genes_diff_score = gene_score_check[gene_score_check > 1]
print(f"Genes where patients have different raw_scores: {len(genes_diff_score)}")

# See examples
if len(genes_diff_score) > 0:
    example_gene = genes_diff_score.index[0]
    print(f"\nExample - {example_gene}:")
    print(all_patient_genes[all_patient_genes['gene_name'] == example_gene][['patient', 'gene_name', 'raw_score']])

genes_with_diff_scores = gene_score_check[gene_score_check > 1].index.tolist()
print(f"Total genes with different raw scores: {len(genes_with_diff_scores)}")
print(genes_with_diff_scores)
