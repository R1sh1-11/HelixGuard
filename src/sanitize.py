import pandas as pd

def sanitize_genome(df: pd.DataFrame, lookup_func) -> pd.DataFrame:
    """
    Iterates through the dataframe, checks each rsID against the blocklist via lookup_func,
    and replaces flagged genotypes with a safe reference value.
    """
    print("[Sanitizer] Starting genomic sanitization process...")
    sanitized_count = 0
    
    # Create a copy so we don't modify the original dataframe unexpectedly
    clean_df = df.copy()
    
    # Loop through every row in our dataset
    for index, row in clean_df.iterrows():
        rsid = row['rsid']
        
        # Call lookup function: expected to return (is_flagged, reference_genotype)
        is_flagged, ref_genotype = lookup_func(rsid)
        
        if is_flagged:
            # Check what the original genotype was before overwriting
            original_genotype = row['genotype']
            
            # If the original genotype already matches the safe reference, no need to overwrite
            if original_genotype != ref_genotype:
                clean_df.at[index, 'genotype'] = ref_genotype
                sanitized_count += 1
                print(f"  [!] Flagged SNP found! Redacted {rsid}: {original_genotype} -> {ref_genotype}")

    print(f"[Sanitizer] Sanitization complete. Redacted {sanitized_count} sensitive SNPs.")
    return clean_df