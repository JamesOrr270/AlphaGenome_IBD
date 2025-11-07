import requests
import pandas as pd

def get_DNA_sequence_UCSC(genome,chrom,start,end):
    URL = f'https://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chrom};start={start};end={end}'
    response = requests.get(URL)
    sequence = response.json()
    return sequence['dna']

def get_CSV_SNP_data(filename):
    SNP_data = pd.read_csv(filename,sep="_",header=None,names=['variant_id','CHROM','POS','REF','ALT']).drop_duplicates().fillna('').replace('-','')
    return (SNP_data)

# def create_SNPeBot_file(transcription_factors,df):


IBD_SNPs = get_CSV_SNP_data('Data/IBD_SNPS.txt')

IBD_SNPs['long_sequence'] = IBD_SNPs.apply(
    lambda row: get_DNA_sequence_UCSC('hg38', row['CHROM'], row['POS']-26, row['POS']+25),
    axis=1
)
IBD_SNPs.to_csv('SNPeBot/IBD_SNPs.csv',index=False)

