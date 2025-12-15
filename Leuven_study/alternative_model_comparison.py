import pandas as pd
import requests
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

tqdm.pandas()

def uniprot_to_gene_name(accession):
    """
    Convert UniProt accession ID to gene name via UniProt REST API.
    
    Args:
        accession (str): UniProt accession ID (e.g., 'P12345')
    
    Returns:
        str: Gene name if found, None otherwise
    """

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

    """
    Convert UniProt accession IDs to gene names for target proteins in the Leuven UC GANDALF dataset.
    
    Reads the affected proteins/TFs/miRs file, converts 'Target' column from UniProt IDs 
    to gene names using the UniProt API, and saves the updated dataframe with a new 
    'Target Gene' column.
    """

    df = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision.txt',
    sep='\t')

    df['Target Gene'] = df['Target'].progress_apply(uniprot_to_gene_name)
    df.to_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name.txt')

def convert_source_data():
    """
    Convert UniProt accession IDs to gene names for source proteins (TFs/miRs) in the dataset.
    
    Reads the file with target gene names already converted, converts 'Source' column 
    from UniProt IDs to gene names using the UniProt API, and saves the dataset with 
    both source and target gene names included.
    """
  
    df = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name.txt',
    sep=',',
    index_col=0)

    df['Source Gene'] = df['Source'].progress_apply(uniprot_to_gene_name)
    df.to_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name_both.txt')

def generate_SNP_impact_lists():

    """
    Generate SNP-to-target gene mappings for each prediction model (FIMO, RSAT, MIRANDA).
    
    Groups target genes by SNP for each model, creating lists of unique genes predicted 
    to be affected by each variant. Falls back to original source identifiers (e.g., miRNA names) 
    when gene names are unavailable.
    
    Returns:
        tuple: Three DataFrames (FIMO, RSAT, MIRANDA), each containing SNPs with their 
               associated lists of target genes
    """

    # Need to make it so that if the source gene has no name to use the orignial code as I think all the time it is just the miRNA
    df = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name_both.txt',
                     sep=',')
    models = ['FIMO','RSAT','MIRANDA']

    df['Source Gene'] = df['Source Gene'].fillna(df['Source'])

    results = [] 


    for model in models:
        model_df = df[df['Interaction_source'] == model]
        model_df = model_df[['Target Gene', 'SNP']]
        
        grouped = model_df.groupby('SNP').agg({
            'Target Gene': lambda x: list(set(x))
        }).reset_index()

        results.append(grouped)
    return results[0], results[1], results[2]
  
def compare_with_AlphaGenome(FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df):
    """
    Compare target gene predictions from FIMO/RSAT/MIRANDA with AlphaGenome predictions across multiple window sizes.
    
    Loads AlphaGenome predictions for four genomic window sizes (16KB, 100KB, 500KB, 1MB), 
    filters for high-confidence predictions (|quantile_score| > 0.99), and merges with 
    traditional model predictions. Calculates overlap metrics and recall scores for each 
    model-window combination.
    
    Args:
        FIMO_grouped_df (DataFrame): FIMO predictions with SNPs and target gene lists
        RSAT_grouped_df (DataFrame): RSAT predictions with SNPs and target gene lists
        MIRANDA_grouped_df (DataFrame): MIRANDA predictions with SNPs and target gene lists
    
    Returns:
        dict: Dictionary with keys like 'FIMO_16KB', 'RSAT_100KB', etc., containing DataFrames 
              with columns: RSID, AlphaGenome (gene list), Model (gene list), Overlap_Count, 
              Model_Gene_Count, AlphaGenome_Gene_count, AG_Recall, Model_Recall
    """

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
    """
    Generate summary statistics for model comparison results across all model-window combinations.
    
    Calculates both per-SNP average recall metrics and overall recall across all genes for each
    model-window pair. AG_Recall represents how many of the traditional model's predicted genes
    were also predicted by AlphaGenome. Model_Recall represents how many of AlphaGenome's 
    predicted genes were captured by the traditional model.
    
    Args:
        merged_results (dict): Dictionary of DataFrames from compare_with_AlphaGenome(), 
                               keyed by 'MODEL_WINDOW' (e.g., 'FIMO_16KB')
    
    Returns:
        DataFrame: Summary table with columns: Model, Window, Total_SNPs, Mean_AG_Recall, 
                   Mean_modal_recall, Overall_AG_Recall, Overall_model_recall
    """
    
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
    """
    Generate SNP-to-source gene mappings for each prediction model (FIMO, RSAT, MIRANDA).
    
    Groups source genes (TFs/miRs) by SNP for each model, creating lists of unique 
    regulatory elements predicted to be affected by each variant. Falls back to original 
    source identifiers (e.g., miRNA names) when gene names are unavailable.
    
    Returns:
        tuple: Three DataFrames (FIMO, RSAT, MIRANDA), each containing SNPs with their 
               associated lists of source genes/regulatory elements
    """
    # Need to make it so that if the source gene has no name to use the orignial code as I think all the time it is just the miRNA
    df = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Data/genetic_data/affected_proteins_TFs_mirs_uc_gandalf_revision_gene_name_both.txt',
                     sep=',')
    models = ['FIMO','RSAT','MIRANDA']

    df['Source Gene'] = df['Source Gene'].fillna(df['Source'])

    results = [] 


    for model in models:
        model_df = df[df['Interaction_source'] == model]
        model_df = model_df[['Source Gene', 'SNP']]
        
        grouped = model_df.groupby('SNP').agg({
            'Source Gene': lambda x: list(set(x))
        }).reset_index()
    
        results.append(grouped)

    return results[0], results[1], results[2]

def compare_source_with_AlphaGenome(FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df):
    """
    Compare source gene (TF/miR) predictions from FIMO/RSAT/MIRANDA with AlphaGenome predictions 
    across multiple window sizes.
    
    Loads AlphaGenome predictions for four genomic window sizes (16KB, 100KB, 500KB, 1MB), 
    filters for high-confidence predictions (|quantile_score| > 0.95), and merges with 
    traditional model source gene predictions. Calculates overlap metrics and AlphaGenome 
    recall scores for each model-window combination.
    
    Args:
        FIMO_grouped_df (DataFrame): FIMO predictions with SNPs and source gene lists
        RSAT_grouped_df (DataFrame): RSAT predictions with SNPs and source gene lists
        MIRANDA_grouped_df (DataFrame): MIRANDA predictions with SNPs and source gene lists
    
    Returns:
        dict: Dictionary with keys like 'FIMO_16KB', 'RSAT_100KB', etc., containing DataFrames 
              with columns: RSID, AlphaGenome (gene list), Model (gene list), Overlap_Count, 
              Model_Gene_Count, AG_Recall
    """

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
        model_name = key.split('_')[0]  
        
        df['Overlap_Count'] = df.apply(
            lambda row: len(set(row['AlphaGenome']) & set(row[model_name])), 
            axis=1
        )
        df['Model_Gene_Count'] = df[model_name].apply(len)
        df['AG_Recall'] = df['Overlap_Count'] / df['Model_Gene_Count']

    return(merged_results)
  
def model_comparison_visualisation(summary_df):

    window_order = ['16KB', '100KB', '500KB', '1MB']
    
    plt.figure(figsize=(10, 6))
    
    for model in ['FIMO', 'RSAT', 'MIRANDA']:
        model_data = summary_df[summary_df['Model'] == model].copy()
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
    
def fishers_exact_test(merged_results):

    results=[]

    for key, df in merged_results.items():
        model_name, window = key.split('_')[0],key.split('_')[1]

        gene_universe_df = pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_{window}_all_scores_RSID.csv')
        unique_gene_count = gene_universe_df['gene_name'].nunique()

        both_predicted = df['Overlap_Count'].sum()
        ag_only = df['AlphaGenome_Gene_count'].sum()-both_predicted
        model_only = df['Model_Gene_Count'].sum()-both_predicted

        total_ag_predictions = df['AlphaGenome_Gene_count'].sum()
        total_model_predictions = df['Model_Gene_Count'].sum()
        
        genes_predicted_by_either = total_ag_predictions + model_only
        neither_predicted = unique_gene_count - genes_predicted_by_either

        table= [[both_predicted,ag_only],
                [model_only,neither_predicted]]

        odds_ratio, p_value = fisher_exact(table, alternative='greater')


        results.append({
            'Model': model_name,
            'Window': window,
            'Both_Predicted': both_predicted,
            'AG_Only': ag_only,
            'Model_Only': model_only,
            'Neither': neither_predicted,
            'Odds_Ratio': odds_ratio,
            'P_value': p_value,
            'Significant': p_value < 0.05
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Leuven_study/Results/model_comparison/fishers_exact_results.csv')
    return(results)

if __name__ == '__main__':
    FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df = generate_SNP_impact_lists()
    gene_comparison = compare_with_AlphaGenome(FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df)
    fishers_exact_results = fishers_exact_test(gene_comparison)


    # summary_df = result_summary(gene_comparison_df)
    # model_comparison_visualisation(summary_df)
