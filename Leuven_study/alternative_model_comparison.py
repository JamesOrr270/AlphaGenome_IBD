import pandas as pd
import requests
from tqdm import tqdm
import matplotlib.pyplot as plt

tqdm.pandas()

def uniprot_to_gene_name(accession):
        try:
        
            url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
            response = requests.get(url)
            
            if response.ok:
                data = response.json()
                gene_name = data['genes'][0]['geneName']['value']
                return gene_name
            return None
        except Exception:
            print(f'Failed: {accession}')

def convert_target_data():
    df = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision.txt',
    sep='\t')

    df['Target Gene'] = df['Target'].progress_apply(uniprot_to_gene_name)
    df.to_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name.txt')

def convert_source_data():
    df = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name.txt',
    sep=',',
    index_col=0)

    df['Source Gene'] = df['Source'].progress_apply(uniprot_to_gene_name)
    df.to_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name_both.txt')

def generate_SNP_impact_lists():
    # Need to make it so that if the source gene has no name to use the orignial code as I think all the time it is just the miRNA
    df = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name_both.txt',
                     sep=',')
    models = ['FIMO','RSAT','MIRANDA']

    df['Source Gene'] = df['Source Gene'].fillna(df['Source'])

    results = []  # Store all three results


    for model in models:
        model_df = df[df['Interaction_source'] == model]
        model_df = model_df[['Target Gene', 'SNP']]
        
        grouped = model_df.groupby('SNP').agg({
            'Target Gene': lambda x: list(set(x))
        }).reset_index()

        results.append(grouped)

    return results[0], results[1], results[2]
  
def compare_with_AlphaGenome(FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df):

    # FIMO only predicts TF binding

    dfs = {
        '16KB': pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_16KB_all_scores_RSID.csv'),
        '100KB': pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_100KB_all_scores_RSID.csv'),
        '500KB': pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_500KB_all_scores_RSID.csv'),
        '1MB': pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_1MB_all_scores_RSID.csv')
    }

    genes = {}
    for window, df in dfs.items():
        df = df[(df['quantile_score'] > 0.99) | (df['quantile_score'] < -0.99)]
        df = df[['gene_name', 'RSID']]
        df = df.groupby('RSID')['gene_name'].apply(lambda x: list(set(x))).reset_index()
        df['RSID'] = df['RSID'].str.upper()
        genes[window] = df

    models = {
        'FIMO': FIMO_grouped_df,
        'RSAT': RSAT_grouped_df,
        'MIRANDA': MIRANDA_grouped_df
    }

    merged_results = {}
    for model_name, model_df in models.items():
        for window, gene_df in genes.items():
            merged = pd.merge(gene_df, model_df, left_on='RSID', right_on='SNP')
            merged = merged[['RSID', 'gene_name', 'Target Gene']]
            merged = merged.rename(columns={'gene_name': 'AlphaGenome', 'Target Gene': model_name})
            merged_results[f'{model_name}_{window}'] = merged

    for key, df in merged_results.items():
        model_name = key.split('_')[0]  # Extract FIMO, RSAT, or MIRANDA
        
        df['Overlap_Count'] = df.apply(
            lambda row: len(set(row['AlphaGenome']) & set(row[model_name])), 
            axis=1
        )
        df['Model_Gene_Count'] = df[model_name].apply(len)
        df['AlphaGenome_Gene_count'] = df['AlphaGenome'].apply(len)
        df['AG_Recall'] = df['Overlap_Count'] / df['Model_Gene_Count']
        df['Model_Recall'] = df['Overlap_Count']/df['AlphaGenome_Gene_count']

    return(merged_results)

def result_summary(merged_results):
    
    summary_data = []
    
    for key, df in merged_results.items():
        model_name, window = key.split('_')[0], key.split('_')[1]
        
        total_snps = len(df)
        mean_AG_recall = df['AG_Recall'].mean()
        mean_modal_recall = df['Model_Recall'].mean()
        
        total_model_genes = df['Model_Gene_Count'].sum()
        total_overlap_genes = df['Overlap_Count'].sum()
        total_AG_genes = df['AlphaGenome_Gene_count'].sum()  

        overall_AG_recall = total_overlap_genes / total_model_genes if total_model_genes > 0 else 0
        overall_Model_recall = total_overlap_genes / total_AG_genes if total_AG_genes > 0 else 0  # Add this

        
        summary_data.append({
            'Model': model_name,
            'Window': window,
            'Total_SNPs': total_snps,
            'Mean_AG_Recall': mean_AG_recall,
            'Mean_modal_recall': mean_modal_recall,
            'Overall_AG_Recall': overall_AG_recall,
            'Overall_model_recall': overall_Model_recall
        })
    
    summary_df = pd.DataFrame(summary_data)
    
    summary_df.to_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Results/model_comparison/summary_table_0.99')
    print(summary_df)
    return(summary_df)
def generate_source_SNP_impact_lists():
    # Need to make it so that if the source gene has no name to use the orignial code as I think all the time it is just the miRNA
    df = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name_both.txt',
                     sep=',')
    models = ['FIMO','RSAT','MIRANDA']

    df['Source Gene'] = df['Source Gene'].fillna(df['Source'])

    results = []  # Store all three results


    for model in models:
        model_df = df[df['Interaction_source'] == model]
        model_df = model_df[['Source Gene', 'SNP']]
        
        grouped = model_df.groupby('SNP').agg({
            'Source Gene': lambda x: list(set(x))
        }).reset_index()

        results.append(grouped)

    return results[0], results[1], results[2]

def compare_source_with_AlphaGenome(FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df):

    # FIMO only predicts TF binding

    dfs = {
        '16KB': pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_16KB_all_scores_RSID.csv'),
        '100KB': pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_100KB_all_scores_RSID.csv'),
        '500KB': pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_500KB_all_scores_RSID.csv'),
        '1MB': pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_1MB_all_scores_RSID.csv')
    }

    genes = {}
    for window, df in dfs.items():
        df = df[(df['quantile_score'] > 0.95) | (df['quantile_score'] < -0.95)]
        df = df[['gene_name', 'RSID']]
        df = df.groupby('RSID')['gene_name'].apply(lambda x: list(set(x))).reset_index()
        df['RSID'] = df['RSID'].str.upper()
        genes[window] = df

    models = {
        'FIMO': FIMO_grouped_df,
        'RSAT': RSAT_grouped_df,
        'MIRANDA': MIRANDA_grouped_df
    }

    merged_results = {}
    for model_name, model_df in models.items():
        for window, gene_df in genes.items():
            merged = pd.merge(gene_df, model_df, left_on='RSID', right_on='SNP')
            merged = merged[['RSID', 'gene_name', 'Source Gene']]
            merged = merged.rename(columns={'gene_name': 'AlphaGenome', 'Source Gene': model_name})
            merged_results[f'{model_name}_{window}'] = merged

    for key, df in merged_results.items():
        model_name = key.split('_')[0]  # Extract FIMO, RSAT, or MIRANDA
        
        df['Overlap_Count'] = df.apply(
            lambda row: len(set(row['AlphaGenome']) & set(row[model_name])), 
            axis=1
        )
        df['Model_Gene_Count'] = df[model_name].apply(len)
        df['AG_Recall'] = df['Overlap_Count'] / df['Model_Gene_Count']

    return(merged_results)
  
def model_comparison_visualisation(summary_df):
    # Define window order
    window_order = ['16KB', '100KB', '500KB', '1MB']
    
    # Create figure
    plt.figure(figsize=(10, 6))
    
    # Plot each model
    for model in ['FIMO', 'RSAT', 'MIRANDA']:
        model_data = summary_df[summary_df['Model'] == model].copy()
        # Sort by window size
        model_data['Window'] = pd.Categorical(model_data['Window'], categories=window_order, ordered=True)
        model_data = model_data.sort_values('Window')
        
        plt.plot(model_data['Window'], model_data['Mean_AG_Recall'], 
                marker='o', linewidth=2, markersize=8, label=model)
    
    plt.xlabel('Window Size', fontsize=12)
    plt.ylabel('Mean Recall', fontsize=12)
    plt.title('AlphaGenome Mean Recall by Model and Window Size', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Results/model_comparison/mean_recall_by_window.png', dpi=300)
    plt.show()
    

if __name__ == '__main__':
    FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df = generate_SNP_impact_lists()
    gene_comparison_df = compare_with_AlphaGenome(FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df)
    summary_df = result_summary(gene_comparison_df)
    model_comparison_visualisation(summary_df)
