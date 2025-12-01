import requests
import time
from pathlib import Path

def get_snp_info(rsid):
    """
    Fetch SNP information from Ensembl REST API
    Returns: tuple of (chromosome, position, ref_allele, list_of_alt_alleles) or None if failed
    """
    # Remove 'rs' prefix if present and clean up
    snp_id = rsid.strip().lower()
    if not snp_id.startswith('rs'):
        snp_id = 'rs' + snp_id
    
    url = f"https://rest.ensembl.org/variation/human/{snp_id}?"
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(url, headers=headers)
        
        if response.status_code == 200:
            data = response.json()
            
            # Get the first mapping (most common)
            if 'mappings' in data and len(data['mappings']) > 0:
                mapping = data['mappings'][0]
                chromosome = mapping.get('seq_region_name', 'NA')
                # Add 'chr' prefix to chromosome
                if chromosome != 'NA':
                    chromosome = f"chr{chromosome}"
                position = mapping.get('start', 'NA')
                
                # Get alleles from the mapping
                allele_string = mapping.get('allele_string', '')
                
                if allele_string and '/' in allele_string:
                    alleles = allele_string.split('/')
                    ref_allele = alleles[0]
                    # Get all alternate alleles as a list
                    alt_alleles = alleles[1:]
                    
                    return (chromosome, position, ref_allele, alt_alleles)
                else:
                    # Fallback: try to get from ancestral_allele and minor_allele
                    ref_allele = data.get('ancestral_allele', 'NA')
                    alt_allele = data.get('minor_allele', 'NA')
                    return (chromosome, position, ref_allele, [alt_allele])
        
        return None
    
    except Exception as e:
        print(f"Error fetching {rsid}: {e}")
        return None

def process_rsid_file(input_file, output_file, delay=0.2):
    """
    Process a file of RSIDs and create output with extended format
    Creates one line per alternate allele
    
    Parameters:
    - input_file: path to input file with RSIDs (one per line)
    - output_file: path to output file
    - delay: delay between API calls in seconds (to avoid rate limiting)
    """
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        success_count = 0
        fail_count = 0
        total_rows = 0
        
        for line_num, line in enumerate(infile, 1):
            rsid = line.strip()
            
            # Skip empty lines or comments
            if not rsid or rsid.startswith('#'):
                continue
            
            print(f"Processing {rsid}... ({line_num})")
            
            # Get SNP information
            snp_info = get_snp_info(rsid)
            
            if snp_info:
                chrom, pos, ref, alt_list = snp_info
                if ref != 'NA' and alt_list and alt_list[0] != 'NA':
                    # Create one line for each alternate allele
                    for alt in alt_list:
                        extended_format = f"{rsid}_{chrom}_{pos}_{ref}_{alt}"
                        outfile.write(f"{extended_format}\n")
                        total_rows += 1
                    
                    if len(alt_list) > 1:
                        print(f"  ✓ {rsid}_{chrom}_{pos}_{ref}_[{'/'.join(alt_list)}] ({len(alt_list)} rows)")
                    else:
                        print(f"  ✓ {rsid}_{chrom}_{pos}_{ref}_{alt_list[0]}")
                    success_count += 1
                else:
                    extended_format = f"{rsid}_{chrom}_{pos}_NA_NA"
                    outfile.write(f"{extended_format}\n")
                    print(f"  ⚠ Position found but alleles missing")
                    fail_count += 1
                    total_rows += 1
            else:
                # Write NA if information couldn't be fetched
                extended_format = f"{rsid}_NA_NA_NA_NA"
                outfile.write(f"{extended_format}\n")
                print(f"  ✗ Failed to fetch information")
                fail_count += 1
                total_rows += 1
            
            # Delay to avoid rate limiting
            time.sleep(delay)
    
    print(f"\n{'='*50}")
    print(f"Processing complete!")
    print(f"Successfully processed: {success_count} SNPs")
    print(f"Failed or incomplete: {fail_count} SNPs")
    print(f"Total rows written: {total_rows}")
    print(f"Results written to {output_file}")
    print(f"{'='*50}")


if __name__ == "__main__":
    patient_list = ['366','481','880','937','955','1782','1914','2376','2634','3146','3365','3670','3771','3792','4133','4572','5513','5517','6030','6684','7051','7148','7194','7645','7689','7748','7951','8193','8573','8660','8691','8842','8864','9165','9442','9608','9971','10097','10485']

    for patient in patient_list:
    
        input_file = f"/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_SNP_lists/{patient}_SNP_list.txt"  # Your input file
        output_file = f"/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/patient_SNP_AG_input/{patient}_AG_SNPs.txt"  # Your output file
        
        # Check if input file exists
        if not Path(input_file).exists():
            print(f"Error: Input file '{input_file}' not found!")
            print("Please create a file with one RSID per line.")
        else:
            process_rsid_file(input_file, output_file)

