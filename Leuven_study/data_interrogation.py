import pandas as pd

patient_list = ['366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']

# patient_metadata = pd.read_excel('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/expression_data/Microarray_Leuven_UC_uninflamed.xlsx')

# print(patient_metadata)

# for patient in patient_list:
#    print(patient_metadata[patient_metadata['VLECCID'] == patient]['Therapy'])

unique_SNPs = set()

for patient in patient_list:
    with open(f'Gene_expression/patient_SNP_lists/{patient}_SNP_list.txt','r') as f:
        for line in f:
            unique_SNPs.add(line.strip())

print(len(unique_SNPs))
print(unique_SNPs)