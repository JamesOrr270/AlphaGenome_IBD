import pandas as pd

SNPs = pd.read_csv('Data/IBD_SNPS.txt',sep='_',header=None, names = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT'])

# Convert to intergers
SNPs['POS'] = SNPs['POS'].astype(int)

# Sort based on chromosome and position
SNPs = SNPs.sort_values(['CHROM','POS']).reset_index(drop=True)

close_SNPs = []

for chrom, chrom_snps in SNPs.groupby('CHROM'):

    for i in range(len(chrom_snps)):
        j = i + 1
        while j < len(chrom_snps):
            distance = chrom_snps.iloc[j]['POS'] - chrom_snps.iloc[i]['POS']
            if distance > 1000000:
                break
            close_SNPs.append({
                'SNP1': chrom_snps.iloc[i]['variant_id'],
                'SNP2': chrom_snps.iloc[j]['variant_id'],
                'Distance': distance,
                'Chromosome':chrom_snps.iloc[j]['CHROM']})
            j+=1

close_SNPs_df = pd.DataFrame(close_SNPs)

close_SNPs_df.to_csv('Results/close_SNPs/close_SNPs_extended.csv', index=False)

