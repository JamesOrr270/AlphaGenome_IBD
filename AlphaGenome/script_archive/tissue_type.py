from alphagenome import colab_utils
from alphagenome.models import dna_client
import pandas as pd

# with open('AlphaGenome/Data/key.txt', 'r') as file:
#     API_key = file.read()

# dna_model = dna_client.create(API_key)

# output_metadata = dna_model.output_metadata(
#     dna_client.Organism.HOMO_SAPIENS
# ).concatenate()

# output_metadata.to_csv('Results/tissue_types/tissue_types.csv',index=False)

SNPs = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/dataset_combination/merged_SNP_dataset.txt',
                   sep='_',
                   names=['RSID','CHROM','POS','REF','ALT','DIS'])

unique_count = SNPs['RSID'].nunique()
print(unique_count)