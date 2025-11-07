# AlphaGenome_IBD

AlphaGenome.py - Takes CSV files containing SNP data and uses AlphaGenome to predict gene expression data

data_prep_pipeline - runs the entire dataprep pipeline

data_sort.py - Takes excel file detailing IBD SNPs and puts them in the correct format for upload to main.py. Also filters out for non-coding SNPs Note: this is tailored to the specific excel file used

tissue_type.py - Lists the tissue types available for use with AlphaGenome

local_SNP.py - Idenifies SNPs that are within 1 MB with each other 

close_SNP_remove.py - Removes SNPs that are close to each other from the main datasets