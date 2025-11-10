import pandas as pd

def read_and_sort(filepath):
    SNPs = pd.read_csv(filepath,sep='_',header=None, names = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT'])
    return(SNPs)

def find_close_SNPs(SNPs):
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
    return(close_SNPs_df)

def write_to_file(close_SNPs_df,filepath):
    close_SNPs_df.to_csv(filepath, index=False)


def main():
    SNPs = read_and_sort('AlphaGenome/Data/IBD_SNPS.txt')

    SNPs['POS'] = SNPs['POS'].astype(int)
    
    close_SNPs_df = find_close_SNPs(SNPs)

    write_to_file(close_SNPs_df,'AlphaGenome/Results/close_SNPs/close_SNPs_extended.csv')

if __name__ == "__main__":
    main()
