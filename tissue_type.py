from alphagenome import colab_utils
from alphagenome.models import dna_client
import pandas as pd

dna_model = dna_client.create('AIzaSyBEl5Qrcby0IUGO2MPUxo62R5y2naLHhDc')

output_metadata = dna_model.output_metadata(
    dna_client.Organism.HOMO_SAPIENS
).concatenate()

output_metadata.to_csv('Results/tissue_types/tissue_types.csv',index=False)