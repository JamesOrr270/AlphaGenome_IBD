# AlphaGenome_IBD

main.py - Takes CSV files containing SNP data and uses AlphaGenome to predict gene expression data

data_sort.py - Takes excel file detailing IBD SNPs and puts them in the correct format for upload to main.py. Note: this is tailored to the specific excel file used

tissue_type.py - Lists the tissue types available for use with AlphaGenome

local_SNP.py - Idenifies SNPs that are within 1 MB with each other 