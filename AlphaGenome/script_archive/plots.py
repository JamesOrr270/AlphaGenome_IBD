import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import ttest_rel, ttest_1samp, shapiro, chi2_contingency, linregress, kstest, anderson

def raw_score_histogram():
    data_16KB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_16KB_all_scores.csv')
    data_100KB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_100KB_all_scores.csv')
    data_500KB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_500KB_all_scores.csv')
    data_1MB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_1MB_all_scores.csv')

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    datasets = [data_16KB, data_100KB, data_500KB, data_1MB]
    titles = ['16KB Window', '100KB Window', '500KB Window', '1MB Window']

    axes_flat = axes.flatten()

    all_data = pd.concat([data['raw_score'] for data in datasets])
    max_abs = max(abs(all_data.min()), abs(all_data.max()))

    for idx, (data, title) in enumerate(zip(datasets, titles)):
        ax = axes_flat[idx]
        
        n, bins, patches = ax.hist(data['raw_score'], bins='auto', edgecolor='black', color='#003E74', alpha=0.7)
        ax.set_xlabel('Raw Score', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        
        ax.set_xlim(-max_abs, max_abs)
        
        ax_inset = ax.inset_axes([0.10, 0.6, 0.35, 0.35])
        
        ax_inset.hist(data['raw_score'], bins=bins, edgecolor='black', color='#003E74', alpha=0.7)
     
        zoom_range = max_abs * 0.03  
        ax_inset.set_xlim(-zoom_range, zoom_range)
        
        y_max = ax.get_ylim()[1]
        ax_inset.set_ylim(0, y_max * 0.3)  
        
        ax_inset.grid(True, alpha=0.3, linestyle='--')
        ax_inset.tick_params(labelsize=8)
        
        ax_inset.patch.set_edgecolor('black')
        ax_inset.patch.set_linewidth(1.5)

    plt.tight_layout()
    plt.show()

def significant_genes_bar_chart():

    data_16KB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_16KB_all_scores.csv')
    data_100KB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_100KB_all_scores.csv')
    data_500KB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_500KB_all_scores.csv')
    data_1MB = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_1MB_all_scores.csv')

    sig_data_16KB = data_16KB[(data_16KB['quantile_score'] > 0.99) | (data_16KB['quantile_score'] < -0.99)]
    sig_data_100KB = data_100KB[(data_100KB['quantile_score'] > 0.99) | (data_100KB['quantile_score'] < -0.99)]
    sig_data_500KB = data_500KB[(data_500KB['quantile_score'] > 0.99) | (data_500KB['quantile_score'] < -0.99)]
    sig_data_1MB = data_1MB[(data_1MB['quantile_score'] > 0.99) | (data_1MB['quantile_score'] < -0.99)]
    
    print(sig_data_16KB)

    datasets = [data_16KB, data_100KB, data_500KB, data_1MB]
    titles = ['16KB Window', '100KB Window', '500KB Window', '1MB Window']

if __name__ == '__main__':
    significant_genes_bar_chart()
