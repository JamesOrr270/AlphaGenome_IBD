"""

Problem:
 - Issues the genes are exactly the same in the matrix, this could be an issue in the data intergration, alternatively this could be an issue in the methodology as taking from a database of SNPs
 so each instance of the SNP will only produce a specified log fold change and will therefore not be different unless there are multiple SNPs targeting the same gene
 - There are some columns that are repeated i.e. A few patients have multiple columns - Have averaged these columns

To do:
- Figure out what to do with the duplicates for the microarray data


"""
import pandas as pd


# Preparing the dataset by restructuing for only the required data, ranaming the columns from the celfile name to the corresponding patient and by meaninng any patient and probe repeats

normalised_expression_data = pd.read_csv('Gene_expression/expression_data/palmieri_annotated_expression.csv')
normalised_expression_data = normalised_expression_data.dropna()
normalised_expression_data = normalised_expression_data.drop(['PROBEID','GENENAME'],axis=1)
normalised_expression_data = normalised_expression_data.set_index('SYMBOL')
patient_metadata = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/expression_data/Microarray_Leuven_UC_uninflamed.txt', sep = '\t')
filename_to_vleccid = dict(zip(patient_metadata['File_name'], patient_metadata['VLECCID']))
normalised_expression_data.rename(columns=filename_to_vleccid, inplace=True)
normalised_expression_data = normalised_expression_data.T.groupby(level=0).mean().T
normalised_expression_data = normalised_expression_data.groupby(level=0).mean()


patient_list = ['366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']
file_sizes = ['16KB','100KB','500KB','1MB']

alphagenome_total_gene_list = []

# This section finds the common genes between the two datasets
for patient in patient_list:
    for size in file_sizes:
        alphaGenome_expression = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_gene_data_AG/{patient}_{size}_gene_list.csv')
        alphagenome_total_gene_list.extend(alphaGenome_expression['gene_name'].tolist())

# This line might not be needed anymore
unique_nomalised_expression_data = set(normalised_expression_data.index)
unique_alphagenome_total_gene_list = list(set(alphagenome_total_gene_list))

common_genes = list(set(unique_alphagenome_total_gene_list) & unique_nomalised_expression_data)

# Makes the dataframe containing the expression data for the common genes between the datasets

expression_matrices = {}

for size in file_sizes:
    expression_matrices[size] = pd.DataFrame(index=common_genes)

    for patient in patient_list:
        # Adding AlphaGenome Data
        alphaGenome_expression = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_gene_data_AG/{patient}_{size}_gene_list.csv')
        alphaGenome_expression = alphaGenome_expression.set_index('gene_name')
        column_name_AG = f'{patient}_AG'
        expression_matrices[size][column_name_AG] = alphaGenome_expression['raw_score'].reindex(common_genes)
        
        # Adding microarray expression data
        column_name_MA = f'{patient}_MA'
        patient_id = int(patient)  # Convert string to integer
        if patient_id in normalised_expression_data.columns:
            expression_matrices[size][column_name_MA] = normalised_expression_data[patient_id].reindex(common_genes)
        else:
            print(f"Patient {patient} not found in expression data columns")
    
    expression_matrices[size].to_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/test/{size}_expression_matrix.csv')