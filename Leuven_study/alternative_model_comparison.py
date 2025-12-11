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
    
    FIMO_df = df[df['Interaction_source']=='FIMO']
    RSAT_df = df[df['Interaction_source']=='RSAT']
    MIRANDA_df = df[df['Interaction_source']=='MIRANDA']

    FIMO_df = FIMO_df[['Source Gene','Target Gene', 'SNP']]
    RSAT_df_df = RSAT_df_df[['Source Gene','Target Gene', 'SNP']]
    MIRANDA_df_df = MIRANDA_df[['Source Gene','Target Gene', 'SNP']]

if __name__ == '__main__':
    convert_source_data()