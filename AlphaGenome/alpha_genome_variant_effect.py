"""
This script uses the AlphaGenome API to make gene expression predictions.

Input - SNP list in the form variantID_CHROM_POS_REF_ALT and name of the output file

Output - two CSV files containing predictions. The first is all predictions for IBD (i.e. in colon) and the secound is Predictions filtered by gene types for protein coding, miRNA and lncRNA

Parameters - Can change the size of input sequence between 16Kb and 1MB
"""



from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
import pandas as pd
from tqdm import tqdm
import numpy as np
import sys


# Intialising the model and also inserting my API key
with open('AlphaGenome/Data/key.txt', 'r') as file:
    API_key = file.read()
  
dna_model = dna_client.create(API_key)

# SNP_data needs to be from HG38 human genome or mouse mm10 genome, I have dropped any variants with none or -, should find a way to not do this
SNP_data = pd.read_csv(sys.argv[1],sep="_",header=None,names=['variant_id','CHROM','POS','REF','ALT','TYPE']).drop_duplicates().fillna('').replace('-','')
SNP_data = SNP_data[['variant_id','CHROM','POS','REF','ALT']]

# Checking weather all of the required columns are presant
required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']

for column in required_columns:
  if column not in SNP_data.columns:
    raise ValueError(f'File is missing required column: {column}.')
 
#  Choose the sequence length. Output the effect of the SNP on all relavant genes within this area. Options are 16KB, 100KB, 500KB, 1MB
sequence_length = '100KB'

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

# Tidy scoers
df_scores = variant_scorers.tidy_scores(results)  

# Filtering by 2SDs from the mean to find significant results
significant_df_scores = df_scores[
  (df_scores['raw_score'] > (df_scores['raw_score'].mean()+ (2*df_scores['raw_score'].std())))|
  (df_scores['raw_score']< df_scores['raw_score'].mean()-(2*df_scores['raw_score'].std()))]

# Filter scores

IBD_specific_significant_scores = significant_df_scores[significant_df_scores['biosample_name'].isin(['colonic mucosa','transverse colon','sigmoid colon','mucosa of descending colon','left colon'])]
gene_type_IBD_specific_significant_scores = IBD_specific_significant_scores[IBD_specific_significant_scores['gene_type'].isin(['protein_coding','miRNA','lncRNA'])]

# Save Files
IBD_specific_significant_scores.to_csv(f'AlphaGenome/Results/AlphaGenome/{sys.argv[2]}_all_scores.csv', index=False)
gene_type_IBD_specific_significant_scores.to_csv(f'AlphaGenome/Results/AlphaGenome/{sys.argv[2]}_filtered_scores.csv', index=False)
