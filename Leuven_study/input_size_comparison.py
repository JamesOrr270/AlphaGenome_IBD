import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


patient_list = ['366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']
file_sizes = ['16KB','100KB','500KB','1MB']

def compare_overlapping_genes_patients():
    total_16KB = pd.DataFrame()
    total_100KB = pd.DataFrame()
    total_500KB = pd.DataFrame()
    total_1MB = pd.DataFrame()
    
    for patient in patient_list:
        df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_alphagenome_results/{patient}_AG_results_16KB_non_significant.csv')
        df = df[(df['quantile_score']>0.99) | (df['quantile_score']<-0.99)]
        df = df['gene_name']
        total_16KB = pd.concat([total_16KB,df])

    for patient in patient_list:
        df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_alphagenome_results/{patient}_AG_results_100KB_non_significant.csv')
        df = df[(df['quantile_score']>0.99) | (df['quantile_score']<-0.99)]
        df = df['gene_name']
        total_100KB = pd.concat([total_100KB,df])
    
    for patient in patient_list:
        df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_alphagenome_results/{patient}_AG_results_500KB_non_significant.csv')
        df = df[(df['quantile_score']>0.99) | (df['quantile_score']<-0.99)]
        df = df['gene_name']
        total_500KB = pd.concat([total_500KB,df])
    
    for patient in patient_list:
        df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_alphagenome_results/{patient}_AG_results_1MB_non_significant.csv')
        df = df[(df['quantile_score']>0.99) | (df['quantile_score']<-0.99)]
        df = df['gene_name']
        total_1MB = pd.concat([total_1MB,df])

    total_16KB = total_16KB.drop_duplicates()
    total_100KB = total_100KB.drop_duplicates()
    total_500KB = total_500KB.drop_duplicates()
    total_1MB = total_1MB.drop_duplicates()

    print(f'No of Genes 16KB: {len(total_16KB)}')
    print(f'No of Genes 100KB: {len(total_100KB)}')
    print(f'No of Genes 500KB: {len(total_500KB)}')
    print(f'No of Genes 1MB: {len(total_1MB)}')

    #  Produces a matrix in which details what percentage of the genes in the x-axis are in the genes in the y-axis i.e. first row is what percentage of genes in 16KB are in the other sizes

    overlap_matrix = pd.DataFrame(columns=file_sizes, index=file_sizes)

    overlap_matrix.loc['16KB', '100KB'] = (total_16KB['gene_name'].isin(total_100KB['gene_name']).sum()) / (len(total_16KB['gene_name'])) * 100
    overlap_matrix.loc['16KB', '500KB'] = (total_16KB['gene_name'].isin(total_500KB['gene_name']).sum()) / (len(total_16KB['gene_name'])) * 100
    overlap_matrix.loc['16KB', '1MB'] = (total_16KB['gene_name'].isin(total_1MB['gene_name']).sum()) / (len(total_16KB['gene_name'])) * 100

    overlap_matrix.loc['100KB', '16KB'] = (total_100KB['gene_name'].isin(total_16KB['gene_name']).sum()) / (len(total_100KB['gene_name'])) * 100
    overlap_matrix.loc['100KB', '500KB'] = (total_100KB['gene_name'].isin(total_500KB['gene_name']).sum()) / (len(total_100KB['gene_name'])) * 100
    overlap_matrix.loc['100KB', '1MB'] = (total_100KB['gene_name'].isin(total_1MB['gene_name']).sum()) / (len(total_100KB['gene_name'])) * 100

    overlap_matrix.loc['500KB', '16KB'] = (total_500KB['gene_name'].isin(total_16KB['gene_name']).sum()) / (len(total_500KB['gene_name'])) * 100
    overlap_matrix.loc['500KB', '100KB'] = (total_500KB['gene_name'].isin(total_100KB['gene_name']).sum()) / (len(total_500KB['gene_name'])) * 100
    overlap_matrix.loc['500KB', '1MB'] = (total_500KB['gene_name'].isin(total_1MB['gene_name']).sum()) / (len(total_500KB['gene_name'])) * 100

    overlap_matrix.loc['1MB', '16KB'] = (total_1MB['gene_name'].isin(total_16KB['gene_name']).sum()) / (len(total_1MB['gene_name'])) * 100
    overlap_matrix.loc['1MB', '100KB'] = (total_1MB['gene_name'].isin(total_100KB['gene_name']).sum()) / (len(total_1MB['gene_name'])) * 100
    overlap_matrix.loc['1MB', '500KB'] = (total_1MB['gene_name'].isin(total_500KB['gene_name']).sum()) / (len(total_1MB['gene_name'])) * 100
   
    print('\n')
    print(overlap_matrix)


def compare_overlapping_genes():
    total_16KB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_16KB_all_scores.csv')
    total_100KB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_100KB_all_scores.csv')
    total_500KB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_500KB_all_scores.csv')
    total_1MB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_1MB_all_scores.csv')

    total_16KB = total_16KB[(total_16KB['quantile_score']>0.99) | (total_16KB['quantile_score']<-0.99)]
    total_100KB = total_100KB[(total_100KB['quantile_score']>0.99) | (total_100KB['quantile_score']<-0.99)]
    total_500KB = total_500KB[(total_500KB['quantile_score']>0.99) | (total_500KB['quantile_score']<-0.99)]
    total_1MB = total_1MB[(total_1MB['quantile_score']>0.99) | (total_1MB['quantile_score']<-0.99)]
    
    total_16KB = total_16KB['gene_name'].drop_duplicates()
    total_100KB = total_100KB['gene_name'].drop_duplicates()
    total_500KB = total_500KB['gene_name'].drop_duplicates()
    total_1MB = total_1MB['gene_name'].drop_duplicates()

    print(f'No of Genes 16KB: {len(total_16KB)}')
    print(f'No of Genes 100KB: {len(total_100KB)}')
    print(f'No of Genes 500KB: {len(total_500KB)}')
    print(f'No of Genes 1MB: {len(total_1MB)}')

    overlap_matrix = pd.DataFrame(columns=file_sizes, index=file_sizes)

    overlap_matrix.loc['16KB', '100KB'] = (total_16KB.isin(total_100KB).sum()) / (len(total_16KB)) * 100
    overlap_matrix.loc['16KB', '500KB'] = (total_16KB.isin(total_500KB).sum()) / (len(total_16KB)) * 100
    overlap_matrix.loc['16KB', '1MB'] = (total_16KB.isin(total_1MB).sum()) / (len(total_16KB)) * 100

    overlap_matrix.loc['100KB', '16KB'] = (total_100KB.isin(total_16KB).sum()) / (len(total_100KB)) * 100
    overlap_matrix.loc['100KB', '500KB'] = (total_100KB.isin(total_500KB).sum()) / (len(total_100KB)) * 100
    overlap_matrix.loc['100KB', '1MB'] = (total_100KB.isin(total_1MB).sum()) / (len(total_100KB)) * 100

    overlap_matrix.loc['500KB', '16KB'] = (total_500KB.isin(total_16KB).sum()) / (len(total_500KB)) * 100
    overlap_matrix.loc['500KB', '100KB'] = (total_500KB.isin(total_100KB).sum()) / (len(total_500KB)) * 100
    overlap_matrix.loc['500KB', '1MB'] = (total_500KB.isin(total_1MB).sum()) / (len(total_500KB)) * 100

    overlap_matrix.loc['1MB', '16KB'] = (total_1MB.isin(total_16KB).sum()) / (len(total_1MB)) * 100
    overlap_matrix.loc['1MB', '100KB'] = (total_1MB.isin(total_100KB).sum()) / (len(total_1MB)) * 100
    overlap_matrix.loc['1MB', '500KB'] = (total_1MB.isin(total_500KB).sum()) / (len(total_1MB)) * 100
    
    print('\n')
    print(overlap_matrix)

    x = file_sizes
    y = [len(total_16KB), len(total_100KB), len(total_500KB), len(total_1MB)]

    plt.figure(figsize=(12, 8))

    colors = ['#003E74']
    bars = plt.bar(x, y, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)

    for i, (bar, value) in enumerate(zip(bars, y)):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{value:,}',  
                ha='center', va='bottom', fontsize=25, fontweight='bold')

    plt.xlabel('Input Sequence Size', fontsize=20, fontweight='bold')
    plt.ylabel('Significantly Differentially Expressed Genes', fontsize=20, fontweight='bold')

    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

   
    plt.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.7)

    
    plt.ylim(0, max(y) * 1.1)

    plt.tight_layout()
    plt.show()

    
    annot_array = overlap_matrix.astype(float).values.copy()
    annot_labels = np.empty_like(annot_array, dtype=object)

    for i in range(annot_array.shape[0]):
        for j in range(annot_array.shape[1]):
            if i == j:  # Diagonal
                annot_labels[i, j] = '-'
            else:
                annot_labels[i, j] = f'{annot_array[i, j]:.1f}%'

    plt.figure(figsize=(14, 6))
    sns.heatmap(overlap_matrix.astype(float), 
                annot=annot_labels,  # Changed from True to annot_labels
                fmt='',   # Changed from '.1f' to empty string
                cmap='YlOrRd',  
                cbar_kws={'label': 'Gene Overlap (%)'},
                linewidths=0.5,
                linecolor='gray',
                vmin=0, vmax=100,
                annot_kws={'fontsize': 20})  # Added this line

    plt.xlabel('', fontsize=12, fontweight='bold')

    plt.gca().xaxis.tick_top()
    plt.gca().xaxis.set_label_position('top')

    plt.xticks(fontsize=20, fontweight='bold')
    plt.yticks(fontsize=20, fontweight='bold', rotation=0)  

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':

    compare_overlapping_genes()

