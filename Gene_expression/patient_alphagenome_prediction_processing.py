import pandas as pd
import requests
import time
import SNP_affected_genes

patients = ['366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']
file_sizes = ['16KB','100KB','500KB','1MB']


def significant_results(filename):

    df = pd.read_csv(filename,usecols=['variant_id','gene_id','gene_name','gene_type','biosample_name','raw_score','quantile_score'])

    df = df[(df['quantile_score'] > 0.99) | (df['quantile_score'] < -0.99)]

    return(df)

def get_rsid_from_api(chrom, pos, ref, alt):
    """
    Get RSID from Ensembl API based on chromosome, position, ref, and alt alleles
    """
    # Remove 'chr' prefix if present
    chrom = chrom.replace('chr', '')
    
    # Ensembl REST API endpoint
    server = "https://rest.ensembl.org"
    ext = f"/overlap/region/human/{chrom}:{pos}-{pos}?feature=variation"
    
    try:
        response = requests.get(server + ext, headers={"Content-Type": "application/json"})
        
        if response.status_code == 200:
            data = response.json()
            
            # Look through variants at this position
            for variant in data:
                if 'id' in variant:
                    # Get alleles - they're already a list
                    alleles = variant.get('alleles', [])
                    
                    # Check if alleles is a list or string and handle accordingly
                    if isinstance(alleles, str):
                        alleles = alleles.split('/')
                    
                    # Check if alleles match
                    if ref in alleles and alt in alleles:
                        return variant['id']
            
            return None
        else:
            return None
            
    except Exception as e:
        print(f"Error fetching RSID for {chrom}:{pos} {ref}>{alt}: {e}")
        return None

def get_RSID(df):
    split_cols = df['variant_id'].str.split(":", expand=True)
    df['CHR'] = split_cols[0] 
    df['POS'] = split_cols[1].astype(int)
    alleles = split_cols[2].str.split(">", expand=True)
    df['REF'] = alleles[0]
    df['ALT'] = alleles[1]  

    df['REF'] = df['REF'].replace('','-').fillna('-')
    df['ALT'] = df['ALT'].replace('','-').fillna('-')

    df['CHR'] = df['CHR'].astype(str)
    
    # Get RSID for each variant using API
    rsids = []
    for idx, row in df.iterrows():
        rsid = get_rsid_from_api(row['CHR'], row['POS'], row['REF'], row['ALT'])
        rsids.append(rsid)
        
        time.sleep(0.1)
        
        # Progress indicator
        if (idx + 1) % 10 == 0:
            print(f"Processed {idx + 1}/{len(df)} variants")
    
    df['RSID'] = rsids
    
    return df

def get_RSID_without_API(df):
    split_cols = df['variant_id'].str.split(":", expand=True)
    df['CHROM'] = split_cols[0] 
    df['POS'] = split_cols[1].astype(int)
    alleles = split_cols[2].str.split(">", expand=True)
    df['REF'] = alleles[0]
    df['ALT'] = alleles[1]  

    df['REF'] = df['REF'].replace('','-').fillna('-')
    df['ALT'] = df['ALT'].replace('','-').fillna('-')

    df['CHROM'] = df['CHROM'].astype(str)

    SNPs = pd.read_csv('/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/dataset_combination/merged_SNP_dataset.txt',
                        sep='_',
                        names=['variant_id','CHROM','POS','REF','ALT','DIS'])
    SNPs = SNPs.replace('','-').fillna('-')
    
    df_scores = pd.merge(
        df,
        SNPs[['variant_id','CHROM','POS','REF','ALT']],
        on=['CHROM', 'POS', 'REF', 'ALT'],
        how='left'
    )

    return(df_scores)
if __name__ == '__main__':
     for size in file_sizes:
          for patient in patients:
            df = get_RSID_without_API(pd.read_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_{size}_all_scores.csv'))
            df.to_csv(f'/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/AlphaGenome/Results/AlphaGenome/All_Nonsig_SNPS_{size}_all_scores_RSID.csv')
            
       
