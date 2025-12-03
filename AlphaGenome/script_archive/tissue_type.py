from alphagenome import colab_utils
from alphagenome.models import dna_client
import pandas as pd

with open('AlphaGenome/Data/key.txt', 'r') as file:
    API_key = file.read()

dna_model = dna_client.create(API_key)

output_metadata = dna_model.output_metadata(
    dna_client.Organism.HOMO_SAPIENS
).concatenate()

output_metadata.to_csv('Results/tissue_types/tissue_types.csv',index=False)