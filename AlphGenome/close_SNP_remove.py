import pandas as pd

def close_SNPs_load(filepath):
    close_SNPs = pd.read_csv(filepath)
    close_SNPs_list = pd.concat([close_SNPs['SNP1'] , close_SNPs['SNP2']]).to_list()
    return(close_SNPs_list)

def load_SNPs(filepath):
    SNP_list = pd.read_csv(filepath, header = None)[0].values.tolist()
    return(SNP_list)

def remove_close_SNPs(disease_SNPs,close_SNP_list):
    SNPs_removed = [snp for snp in disease_SNPs if not any(substring in str(snp) for substring in close_SNP_list)]
    return(SNPs_removed)

def write_to_file(data,filepath):
    with open(filepath, 'w') as f:
        for item in data:
            f.write(f"{item}\n")

def main():

    close_SNPs = close_SNPs_load('Results/close_SNPs/close_SNPs_extended.csv')

    diseases = ['IBD', 'CD', 'UC']
     
    for disease in diseases:
        
        snps = load_SNPs(f'Data/{disease}_SNPS.txt')
        filtered = remove_close_SNPs(snps, close_SNPs)
        write_to_file(filtered, f'Results/close_SNPs_Removed/{disease}_SNPs_removed_test3.txt') 

if __name__ == '__main__':
    main()
