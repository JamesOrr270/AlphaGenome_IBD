"""
Unit tests for SNP dataset merging functionality.
Tests the merging logic described in the main program comments.
"""
import pandas as pd
import numpy as np
import unittest
import tempfile
import os


class TestSNPMerging(unittest.TestCase):
    
    def setUp(self):
        """Create temporary test files before each test"""
        self.temp_dir = tempfile.mkdtemp()
        self.website_file = os.path.join(self.temp_dir, 'website_SNP_list.txt')
        self.raw_file = os.path.join(self.temp_dir, 'raw_SNP_list.txt')
        self.output_file = os.path.join(self.temp_dir, 'merged_SNP_dataset.txt')
    
    def tearDown(self):
        """Clean up temporary files after each test"""
        for file in [self.website_file, self.raw_file, self.output_file]:
            if os.path.exists(file):
                os.remove(file)
        os.rmdir(self.temp_dir)
    
    def create_test_data(self, raw_data, website_data):
        """Helper function to create test input files"""
        # Write raw SNP data
        with open(self.raw_file, 'w') as f:
            for row in raw_data:
                f.write('_'.join(str(x) if pd.notna(x) else '' for x in row) + '\n')
        
        # Write website SNP data
        with open(self.website_file, 'w') as f:
            for row in website_data:
                f.write('_'.join(str(x) if pd.notna(x) else '' for x in row) + '\n')
    
    def run_merge_process(self):
        """Run the main merging logic"""
        # Load the data
        website_SNPs = pd.read_csv(self.website_file, sep='_', header=None, 
                                   names=['RSID','CHR', 'POS','REF','ALT','TYPE'])
        
        Raw_SNPs = pd.read_csv(self.raw_file, sep='_', header=None, 
                              names=['RSID','CHR','POS','REF','ALT','TYPE'])
        
        # Refine website SNPs
        website_SNP_refined = website_SNPs[['RSID', 'TYPE']].drop_duplicates()
        
        website_SNP_refined = website_SNP_refined.groupby('RSID')['TYPE'].apply(
            lambda x: ' '.join(x.dropna().unique())
        ).reset_index()
        
        website_SNP_refined['TYPE'] = website_SNP_refined['TYPE'].replace('', np.nan)
        
        # Merge datasets
        Raw_SNPs = Raw_SNPs.drop_duplicates()
        merged_SNPs = pd.merge(Raw_SNPs, website_SNP_refined, on='RSID', 
                              how='left', suffixes=('_raw','_web'))
        
        def combine_types(row):
            raw_val = row['TYPE_raw']
            web_val = row['TYPE_web']
            
            if pd.notna(raw_val) and pd.notna(web_val):
                total_types = ' '.join(set(raw_val.split(' ') + web_val.split(' ')))
                return total_types
            elif pd.isna(raw_val) and pd.notna(web_val):
                return web_val
            elif pd.notna(raw_val) and pd.isna(web_val):
                return raw_val
            else:
                return np.nan
        
        merged_SNPs['TYPE'] = merged_SNPs.apply(combine_types, axis=1)
        merged_SNPs = merged_SNPs.drop(columns=['TYPE_raw', 'TYPE_web'])
        merged_SNPs = merged_SNPs.drop_duplicates()
        
        # Write output
        merged_SNPs.to_csv(self.output_file, sep='_', index=None, header=None)
        
        return merged_SNPs
    
    def test_basic_merge_with_type_addition(self):
        """Test that website type data is added to raw SNPs based on RSID"""
        raw_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', np.nan],
            ['rs456', 'CHR2', '2000', 'C', 'T', np.nan]
        ]
        website_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'CD'],
            ['rs456', 'CHR2', '2000', 'C', 'T', 'UC']
        ]
        
        self.create_test_data(raw_data, website_data)
        result = self.run_merge_process()
        
        # Check that types were added
        self.assertEqual(result[result['RSID'] == 'rs123']['TYPE'].values[0], 'CD')
        self.assertEqual(result[result['RSID'] == 'rs456']['TYPE'].values[0], 'UC')
    
    def test_snps_without_type_remain_blank(self):
        """Test that SNPs with no disease type from both files are left blank (NaN)"""
        raw_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', np.nan],
        ]
        website_data = [
            ['rs456', 'CHR2', '2000', 'C', 'T', 'CD'],  # Different RSID
        ]
        
        self.create_test_data(raw_data, website_data)
        result = self.run_merge_process()
        
        # rs123 should have NaN type since it's not in website data
        self.assertTrue(pd.isna(result[result['RSID'] == 'rs123']['TYPE'].values[0]))
    
    def test_multiple_phenotypes_combined(self):
        """Test that SNPs with multiple phenotypes have all phenotypes in TYPE column"""
        raw_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'CD'],
        ]
        website_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'UC'],
            ['rs123', 'CHR1', '1000', 'A', 'G', 'IBD']
        ]
        
        self.create_test_data(raw_data, website_data)
        result = self.run_merge_process()
        
        type_value = result[result['RSID'] == 'rs123']['TYPE'].values[0]
        # Check that all three phenotypes are present
        self.assertIn('CD', type_value)
        self.assertIn('UC', type_value)
        self.assertIn('IBD', type_value)
    
    def test_unique_phenotypes_no_duplicates(self):
        """Test that duplicate phenotypes are removed"""
        raw_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'CD'],
        ]
        website_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'CD'],  # Duplicate CD
            ['rs123', 'CHR1', '1000', 'A', 'G', 'UC']
        ]
        
        self.create_test_data(raw_data, website_data)
        result = self.run_merge_process()
        
        type_value = result[result['RSID'] == 'rs123']['TYPE'].values[0]
        phenotypes = type_value.split(' ')
        
        # Check CD appears only once
        self.assertEqual(phenotypes.count('CD'), 1)
        self.assertIn('UC', phenotypes)
    
    def test_all_raw_snps_preserved(self):
        """Test that all RSIDs from raw file are in the output"""
        raw_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'CD'],
            ['rs456', 'CHR2', '2000', 'C', 'T', np.nan],
            ['rs789', 'CHR3', '3000', 'G', 'A', 'UC']
        ]
        website_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'IBD'],
        ]
        
        self.create_test_data(raw_data, website_data)
        result = self.run_merge_process()
        
        # Check all RSIDs from raw data are present
        self.assertEqual(len(result), 3)
        self.assertIn('rs123', result['RSID'].values)
        self.assertIn('rs456', result['RSID'].values)
        self.assertIn('rs789', result['RSID'].values)
    
    def test_rsid_only_matching(self):
        """Test that matching is based only on RSID, not nucleotide change"""
        raw_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', np.nan],  # A->G
        ]
        website_data = [
            ['rs123', 'CHR1', '1000', 'C', 'T', 'CD'],  # C->T (different nucleotides)
        ]
        
        self.create_test_data(raw_data, website_data)
        result = self.run_merge_process()
        
        # Type should be added despite different nucleotides
        self.assertEqual(result[result['RSID'] == 'rs123']['TYPE'].values[0], 'CD')
    
    def test_raw_type_preserved_when_no_website_match(self):
        """Test that raw SNP type is preserved when there's no website match"""
        raw_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'CD'],
        ]
        website_data = [
            ['rs999', 'CHR9', '9000', 'T', 'A', 'UC'],  # Different RSID
        ]
        
        self.create_test_data(raw_data, website_data)
        result = self.run_merge_process()
        
        # Original CD type should be preserved
        self.assertEqual(result[result['RSID'] == 'rs123']['TYPE'].values[0], 'CD')
    
    def test_duplicate_removal(self):
        """Test that duplicate SNPs are removed"""
        raw_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'CD'],
            ['rs123', 'CHR1', '1000', 'A', 'G', 'CD'],  # Duplicate
        ]
        website_data = [
            ['rs456', 'CHR2', '2000', 'C', 'T', 'UC'],
        ]
        
        self.create_test_data(raw_data, website_data)
        result = self.run_merge_process()
        
        # Should only have one rs123 entry
        self.assertEqual(len(result[result['RSID'] == 'rs123']), 1)
    
    def test_output_file_created(self):
        """Test that the output file is created successfully"""
        raw_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'CD'],
        ]
        website_data = [
            ['rs123', 'CHR1', '1000', 'A', 'G', 'UC'],
        ]
        
        self.create_test_data(raw_data, website_data)
        self.run_merge_process()
        
        # Check output file exists
        self.assertTrue(os.path.exists(self.output_file))
        
        # Check it can be read back
        output_df = pd.read_csv(self.output_file, sep='_', header=None,
                               names=['RSID','CHR','POS','REF','ALT','TYPE'])
        self.assertEqual(len(output_df), 1)


if __name__ == '__main__':
    unittest.main()