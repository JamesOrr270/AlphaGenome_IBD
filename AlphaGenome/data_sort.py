import pandas as pd

# Function that takes the data and puts into into the relavant files for analysis
def load_data(SNP_path,non_coding_types_path):
    SNP_file = pd.read_excel(SNP_path,skiprows=1)
    with  open(non_coding_types_path, 'r') as f:
        non_coding_types = f.read().splitlines()
    return SNP_file, non_coding_types

# Splitting into non-coding SNPs and coding SNPs
def is_non_coding(consequence,non_coding_types):
    if pd.isna(consequence):
        return False
    consequence_list = [item.strip() for item in str(consequence).split(',')]
    return any(item in non_coding_types for item in consequence_list)

# Creates a mask that can filter if coding or non-coding
def filter_non_coding_snps(SNP_file, non_coding_types):
    non_coding_mask = SNP_file['Consequence'].apply(
        lambda x: is_non_coding(x, non_coding_types)
    )
    return SNP_file[non_coding_mask], SNP_file[~non_coding_mask]

# Puts the strings into the correct format
def create_snp_string(row):
    return f"{row.RS_ID}_chr{row.Chr_index}_{row.Pos_index}_{row.A1_index}_{row.A2_index}"

# Splits SNPs into the differnet subtypes
def IBD_type_split(non_coding_SNPs,subtype):
    filtered = non_coding_SNPs[non_coding_SNPs['Phenotype_loci'] == subtype]
    return [create_snp_string(row) for row in filtered.itertuples()]

# Aggregates the lists of different SNPs
def aggregate_snps(non_coding_SNPs):
   
    cd_snps = IBD_type_split(non_coding_SNPs, 'CD')
    uc_snps = IBD_type_split(non_coding_SNPs, 'UC')
    ibd_snps = IBD_type_split(non_coding_SNPs, 'IBD')
    
    all_snps = cd_snps + uc_snps + ibd_snps
    cd_snps_combined = cd_snps + ibd_snps
    uc_snps_combined = uc_snps + ibd_snps
    
 
    return {
        'All': list(dict.fromkeys(all_snps)),
        'CD': list(dict.fromkeys(cd_snps_combined)),
        'UC': list(dict.fromkeys(uc_snps_combined))
    }

# Saves the SNP lists to a file
def write_snp_list(filepath, snp_list):
    with open(filepath, 'w') as file:
        file.write('\n'.join(snp_list) + '\n')

def main(file_SNP,file_non_coding):
    SNP_df, non_coding_types = load_data(file_SNP,file_non_coding)
    
    non_coding_SNPs, coding_SNPs = filter_non_coding_snps(SNP_df, non_coding_types)
    
    snp_collections = aggregate_snps(non_coding_SNPs)
    
    write_snp_list('AlphaGenome/Data/IBD_SNPS.txt', snp_collections['All'])
    write_snp_list('AlphaGenome/Data/CD_SNPS.txt', snp_collections['CD'])
    write_snp_list('AlphaGenome/Data/UC_SNPS.txt', snp_collections['UC'])
    

if __name__ == "__main__":
    
    SNP_file = 'Data/Liu. et al.2023 studys 320 SNPs.xlsx'
    non_coding_types_file = 'Data/non_coding_types.txt'

    main(SNP_file,non_coding_types_file)





