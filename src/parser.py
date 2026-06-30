import os
from snps import SNPs
import pandas as pd

def parse_genome(file_path: str) -> pd.DataFrame:
    #loads raw genome file using the snps library
    #returns pandas DataFrame with rsid, chromosome, position, genotype columns
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Genome file not found at: {file_path}")
        
    print(f"[Parser] Loading genome file: {file_path}...")
    
    #the snps library automatically handles headers and different vendor formats
    parsed_data = SNPs(file_path)
    
    #convert the internal SNPs dictionary/dataframe to a clean DataFrame for processing
    df = parsed_data.snps.reset_index()
    
    #ensure columns match our expected pipeline format
    df.columns = ['rsid', 'chrom', 'pos', 'genotype']
    
    print(f"[Parser] Successfully loaded {len(df)} SNPs.")
    return df