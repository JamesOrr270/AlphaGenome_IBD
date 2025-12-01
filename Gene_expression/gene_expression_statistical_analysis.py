"""
This script analyses the patient specific gene expression data. Due to the limitations in the data does a binary comparison in which it
compares up and downregulation. Takes the output from patient gene_expression_matrix.py as input and outputs results 
into the statistical_test_results folder. 

Functions:
one_tail_analysis() - Performs one tail t-test on results followed by benjamini hochberg multiple testing correction
two_tail_analysis() - Performs two tail t-test on results followed by benjamini hochberg multiple testing correction

Parameters:
- Need to choose between one and two tail test"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import spearmanr, pearsonr, norm
import matplotlib.pyplot as plt
import seaborn as sn
from statsmodels.stats.multitest import multipletests

patient_list = ['366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']
file_sizes = ['16KB','100KB','500KB','1MB']
# Making a matrix with 0s and 1s could be good for t-test

def one_tail_analysis():
    for size in file_sizes:
        df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/gene_expression_matrix_frompatient/{size}_expression_matrix.csv',index_col=0)
        MA_matrix = df[[col for col in df.columns if col.endswith('_MA')]].copy()
        AG_matrix = df[[col for col in df.columns if col.endswith('_AG')]].copy()
        AG_matrix = np.sign(AG_matrix.fillna(0))
        total_matrix = pd.merge(AG_matrix, MA_matrix, left_index=True, right_index=True, how='outer')

        significant_count = 0
        total_tested = 0
        significant_genes = []
        p_values = []
        gene_names = []
        directions = []
        
        for gene in total_matrix.index:
            group_1 = [] #Group 1 is the group with all the non-differentially expressed
            group_2 = [] #Group 2 is the group with differentially expressed
            direction = None

            for patient in patient_list:
                if total_matrix.loc[gene,f'{patient}_AG'] == 0:
                    group_1.append(total_matrix.loc[gene,f'{patient}_MA'])
                else:
                    group_2.append(total_matrix.loc[gene,f'{patient}_MA'])
                    if total_matrix.loc[gene,f'{patient}_AG'] > 0:
                        direction = 'upregulation'
                    elif total_matrix.loc[gene,f'{patient}_AG'] < 0:
                        direction = 'downregulation'

            if len(group_1) < 3 or len(group_2) < 3:
                print(f'For {gene} there are not enough values in either group 1 or group 2')
                continue

            if direction == 'upregulation':
                    t_stat, p_value = stats.ttest_ind(group_1, group_2, alternative='less')
                    total_tested += 1
                    if p_value <0.05:
                        significant_count +=1
                        significant_genes.append(gene)
            elif direction == 'downregulation':
                    t_stat, p_value = stats.ttest_ind(group_1, group_2, alternative='greater')
                    total_tested += 1
                    if p_value < 0.05:
                        significant_count +=1
                        significant_genes.append(gene)
                
            elif direction == None:
                print(f'{gene} are all 0 values')
                continue
            
            p_values.append(p_value)
            gene_names.append(gene)
            directions.append(direction)
        
        
            print(f'--------------{gene}----------------')
            print(f'Testing if group 2 is {direction} compared to group 1')
            print(f'Group 1 Mean:{np.mean(group_1)}')
            print(f'Group 2 Mean:{np.mean(group_2)}')
            print(f"T-statistic: {t_stat}")
            print(f"P-value: {p_value}")


        print(f'Total tested: {total_tested}')
        print(f'Significant results {significant_count}')
        print(f'Significant genes: {significant_genes}')


        # Apply multiple testing correction
        reject, p_corrected, alphacSidak, alphacBonf = multipletests(
            p_values, 
            alpha=0.05, 
            method='fdr_bh'  # Benjamini-Hochberg FDR correction
        )

        # Create results dataframe
        results_df = pd.DataFrame({
            'gene': gene_names,
            'direction': directions,
            'p_value': p_values,
            'p_value_corrected': p_corrected,
            'significant': reject
        })

        results_df.to_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/statistical_test_results/{size}_one_tail.csv')

        # Count significant results
        print(f"Significant genes (uncorrected p < 0.05): {(results_df['p_value'] < 0.05).sum()}")
        print(f"Significant genes (FDR < 0.05): {results_df['significant'].sum()}")

def two_tail_analysis():
      for size in file_sizes:
        df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/gene_expression_matrix_frompatient/{size}_expression_matrix.csv',index_col=0)
        MA_matrix = df[[col for col in df.columns if col.endswith('_MA')]].copy()
        AG_matrix = df[[col for col in df.columns if col.endswith('_AG')]].copy()
        AG_matrix = np.sign(AG_matrix.fillna(0))
        total_matrix = pd.merge(AG_matrix, MA_matrix, left_index=True, right_index=True, how='outer')

        significant_count = 0
        total_tested = 0
        significant_genes = []
        p_values = []
        gene_names = []
        
        for gene in total_matrix.index:
            group_1 = [] #Group 1 is the group with all the non-differentially expressed
            group_2 = [] #Group 2 is the group with differentially expressed

            for patient in patient_list:
                if total_matrix.loc[gene,f'{patient}_AG'] == 0:
                    group_1.append(total_matrix.loc[gene,f'{patient}_MA'])
                else:
                    group_2.append(total_matrix.loc[gene,f'{patient}_MA'])
                    

            if len(group_1) < 3 or len(group_2) < 3:
                print(f'For {gene} there are not enough values in either group 1 or group 2')
                continue

    
            t_stat, p_value = stats.ttest_ind(group_1, group_2)
            total_tested += 1
            if p_value <0.05:
                significant_count +=1
                significant_genes.append(gene)
                
        
            
            p_values.append(p_value)
            gene_names.append(gene)
        
        
        
            print(f'--------------{gene}----------------')
            print(f'Group 1 Mean:{np.mean(group_1)}')
            print(f'Group 2 Mean:{np.mean(group_2)}')
            print(f"T-statistic: {t_stat}")
            print(f"P-value: {p_value}")


        print(f'Total tested: {total_tested}')
        print(f'Significant results {significant_count}')
        print(f'Significant genes: {significant_genes}')


        # Apply multiple testing correction
        reject, p_corrected, alphacSidak, alphacBonf = multipletests(
            p_values, 
            alpha=0.05, 
            method='fdr_bh'  # Benjamini-Hochberg FDR correction
        )

        # Create results dataframe
        results_df = pd.DataFrame({
            'gene': gene_names,
            'p_value': p_values,
            'p_value_corrected': p_corrected,
            'significant': reject
        })

        results_df.to_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/statistical_test_results/{size}_two_tail.csv')

        # Count significant results
        print(f"Significant genes (uncorrected p < 0.05): {(results_df['p_value'] < 0.05).sum()}")
        print(f"Significant genes (FDR < 0.05): {results_df['significant'].sum()}")

if __name__ == '__main__':
    two_tail_analysis()

