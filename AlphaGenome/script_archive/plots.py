import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import ttest_rel, ttest_1samp, shapiro, chi2_contingency, linregress, kstest, anderson



data = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_1MB_all_scores.csv')
# print(data)
# data['raw_score'].hist(bins=30)
# plt.xlabel('Raw Score')
# plt.ylabel('Frequency')
# plt.title('Frequency of Raw Scores')
# plt.show()

# # For discrete/limited values, a bar plot might be better
# plt.figure(figsize=(10, 6))
# data['raw_score'].value_counts().sort_index().plot(kind='bar', edgecolor='black')
# plt.xlabel('Quantile Score')
# plt.ylabel('Count')
# plt.title('Frequency of Quantile Scores')
# plt.show()

# Use more appropriate binning
plt.figure(figsize=(10, 6))
data['raw_score'].hist(bins='auto', edgecolor='#003E74')
plt.xlabel('Raw Score')
plt.ylabel('Frequency')
plt.title('Distribution of Raw Scores')
plt.ylim(0, 50) 
plt.show()

# print(data['raw_score'].describe())
# data['quantile_score'].hist(bins=10)
# plt.xlabel('Quantile Score')
# plt.ylabel('Frequency')
# plt.title('Frequency of Quantile Scores')
# plt.show()

# print(data['quantile_score'].describe())

# data['quantile_score'] = data['quantile_score'].abs()

# data['quantile_score'].hist(bins=40)
# plt.xlabel('Quantile Score')
# plt.ylabel('Frequency')
# plt.title('Frequency of Quantile Scores')
# plt.show()

# print(data['quantile_score'].describe())
