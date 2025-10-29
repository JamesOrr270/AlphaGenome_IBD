import pandas as pd

close_SNPs = pd.read_csv('Results/close_SNPs_extended.csv')
IBD_SNPs = pd.read_csv('Data/IBD_SNPS.txt', header = None)[0].values.tolist()
CD_SNPs = pd.read_csv('Data/CD_SNPS.txt',header = None)[0].values.tolist()
UC_SNPs = pd.read_csv('Data/UC_SNPS.txt',header = None)[0].values.tolist()

close_SNPs_list = pd.concat([close_SNPs['SNP1'] , close_SNPs['SNP2']]).to_list()

IBD_SNPs = [snp for snp in IBD_SNPs if not any(substring in str(snp) for substring in close_SNPs_list)]
CD_SNPs = [snp for snp in CD_SNPs if not any(substring in str(snp) for substring in close_SNPs_list)]
UC_SNPs = [snp for snp in UC_SNPs if not any(substring in str(snp) for substring in close_SNPs_list)]

with open('Results/close_SNPs_Removed/IBD_SNPs_removed.txt', 'w') as f:
    for item in IBD_SNPs:
        f.write(f"{item}\n")

with open('Results/close_SNPs_Removed/CD_SNPs_removed.txt', 'w') as f:
    for item in CD_SNPs:
        f.write(f"{item}\n")

with open('Results/close_SNPs_Removed/UC_SNPs_removed.txt', 'w') as f:
    for item in UC_SNPs:
        f.write(f"{item}\n")