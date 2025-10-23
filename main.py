from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers

import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import numpy as np

# Intialising the model and also inserting my API key
dna_model = dna_client.create('AIzaSyBEl5Qrcby0IUGO2MPUxo62R5y2naLHhDc')

# SNP_data needs to be from HG38 human genome or mouse mm10 genome
SNP_data = pd.read_csv('SNPs_loc/test.txt',sep="_",header=None,names=['variant_id','CHROM','POS','REF','ALT'])
SNP_data = SNP_data.replace(['None','-'],np.nan)
SNP_data = SNP_data.dropna()

# Checking weather all of the required columns are presant
required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']

for column in required_columns:
  if column not in SNP_data.columns:
    raise ValueError(f'VCF file is missing required column: {column}.')
 
#  Need to look at what the best value for this is, should be based on what I want to find out, I chose this one as it is the most comprehensive
sequence_length = '2KB'
# Makes the sequence length in the form that alphaGenome requires
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f'SEQUENCE_LENGTH_{sequence_length}']

# Choosing which tests to predict
score_rna_seq = True
score_cage = False
score_procap = False
score_dnase = False
score_chip_tf = False
score_atac = False
score_chip_histone = False
score_polyadenylation = False
score_splice_sites = False
score_splice_site_usage = False
score_splice_junctions = False

# This determines weather the results are saved 
download_predictions = True

# This puts the species in a form that alphagenome will understand
organism = 'human'
organism_map = {
    'human': dna_client.Organism.HOMO_SAPIENS,
    'mouse': dna_client.Organism.MUS_MUSCULUS,
}
organism = organism_map[organism]

# This section chooses your supported scorer and puts it into a form that AlphaGenome will recognise
scorer_selections = {
    'rna_seq': score_rna_seq,
    'cage': score_cage,
    'procap': score_procap,
    'atac': score_atac,
    'dnase': score_dnase,
    'chip_histone': score_chip_histone,
    'chip_tf': score_chip_tf,
    'polyadenylation': score_polyadenylation,
    'splice_sites': score_splice_sites,
    'splice_site_usage': score_splice_site_usage,
    'splice_junctions': score_splice_junctions,
}

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
selected_scorers = [
    all_scorers[key]
    for key in all_scorers
    if scorer_selections.get(key.lower(), False)
]

# Removes the scorers that are not supported for your organism
unsupported_scorers = [
    scorer
    for scorer in selected_scorers
    if (
        organism.value
        not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
    )
    | (
        (scorer.requested_output == dna_client.OutputType.PROCAP)
        & (organism == dna_client.Organism.MUS_MUSCULUS)
    )
]
if len(unsupported_scorers) > 0:
  print(
      f'Excluding {unsupported_scorers} scorers as they are not supported for'
      f' {organism}.'
  )
  for unsupported_scorer in unsupported_scorers:
    selected_scorers.remove(unsupported_scorer)

# Scores the SNPs
results = []

for i, SNP_row in tqdm(SNP_data.iterrows(), total=len(SNP_data)):
  variant = genome.Variant(
      chromosome=str(SNP_row.CHROM),
      position=int(SNP_row.POS),
      reference_bases=SNP_row.REF,
      alternate_bases=SNP_row.ALT,
      name=SNP_row.variant_id,
  )
  interval = variant.reference_interval.resize(sequence_length)

  variant_scores = dna_model.score_variant(
      interval=interval,
      variant=variant,
      variant_scorers=selected_scorers,
      organism=organism,
  )
  results.append(variant_scores)

# Tidy and filter the scores
df_scores = variant_scorers.tidy_scores(results)
filtered_df_scores = df_scores[df_scores['biosample_name'].isin(['large intestine','small intestine'])]

if download_predictions:
  filtered_df_scores.to_csv('Results/results_test.csv', index=False)
  
