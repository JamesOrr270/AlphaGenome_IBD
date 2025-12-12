import pandas as pd
import requests
from tqdm import tqdm

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
            merged = merged[['RSID', 'gene_name', 'Target Gene']]
            merged = merged.rename(columns={'gene_name': 'AlphaGenome', 'Target Gene': model_name})
            merged_results[f'{model_name}_{window}'] = merged

    print(merged_results['RSAT_16KB'])


if __name__ == '__main__':
    FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df = generate_SNP_impact_lists()
    gene_comparison_df = compare_with_AlphaGenome(FIMO_grouped_df,RSAT_grouped_df,MIRANDA_grouped_df)