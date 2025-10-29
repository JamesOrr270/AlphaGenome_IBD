import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import ttest_rel, ttest_1samp, shapiro, chi2_contingency, linregress, kstest, anderson



data = pd.read_csv('Results/AlphaGenome/IBD_results_29Oct25.csv')

data['raw_score'].hist(bins=30)
plt.xlabel('Raw Score')
plt.ylabel('Frequency')
plt.title('Frequency of Raw Scores')
plt.show()

print(data['raw_score'].describe())
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
