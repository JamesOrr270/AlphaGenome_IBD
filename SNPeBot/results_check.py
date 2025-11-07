import pandas as pd

SNPs = pd.read_csv('SNPeBot/IBD_SNPs.csv')

# Check each sequence length
lengths = SNPs['long_sequence'].str.len()

print(f"Min length: {lengths.min()}")
print(f"Max length: {lengths.max()}")
print(f"All are 51 bp: {(lengths == 51).all()}")

# Show counts of each length
print(lengths.value_counts())

SNPs['base_at_pos_26'] = SNPs['long_sequence'].str[25]

# Check if it matches REF allele
SNPs['matches_ref'] = SNPs['base_at_pos_26'].str.upper() == SNPs['REF'].str.upper()

SNPs.to_csv('SNPeBot/check.csv')

# Summary
print(f"Total SNPs: {len(SNPs)}")
print(f"SNPs where position 26 matches REF: {SNPs['matches_ref'].sum()}")
print(f"SNPs where position 26 does NOT match REF: {(~SNPs['matches_ref']).sum()}")

# Show examples of mismatches
mismatches = SNPs[~SNPs['matches_ref']]
if len(mismatches) > 0:
    print("\nExamples of mismatches:")
    print(mismatches[['CHROM', 'POS', 'REF', 'base_at_pos_26', 'long_sequence']].head())
else:
    print("\nAll SNPs have REF allele at position 26! ✓")