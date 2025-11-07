import data_sort
import local_SNP
import close_SNP_remove

# Put in input file in first and then the non_coding_types in the secound
data_sort.main('Data/Liu. et al.2023 studys 320 SNPs.xlsx','Data/non_coding_types.txt')

local_SNP.main()

close_SNP_remove.main()