import pandas as pd 

patient_list = ['366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']
file_sizes = ['16KB','100KB','500KB','1MB']
for size in file_sizes:

    for patient in patient_list:
        all_gene_expression = pd.DataFrame()
        patient_df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_gene_data_AG/{patient}_{size}_gene_list.csv')
        all_gene_expression = pd.concat([all_gene_expression,patient_df])
        
    print(f'Top 10 for {size}: {all_gene_expression}')


