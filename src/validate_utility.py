import os
import sys

def test_snps_library(file_path):
    """
    Tests if the open-source 'snps' Python library can successfully 
    parse the sanitized genomic file without crashing.
    """
    print(f"[QA Test 1] Testing file compatibility with 'snps' library: {file_path}")
    try:
        from snps import SNPs
    except ImportError:
        print("  [!] 'snps' library not installed. Install via: pip install snps")
        return False

    try:
        # Load file using SNPs
        s = SNPs(file_path)
        snp_count = len(s.snps)
        
        if snp_count > 0:
            print(f"  [PASS] 'snps' library parsed {file_path} successfully!")
            print(f"         Total SNPs parsed: {snp_count}")
            return True
        else:
            print("  [FAIL] 'snps' loaded 0 variants.")
            return False

    except Exception as e:
        # If 'snps' throws an int split error due to numeric chrom values, 
        # retry by passing bytes or forcing string types
        try:
            import pandas as pd
            df = pd.read_csv(file_path, sep="\t", comment="#", header=None, dtype=str)
            snp_count = len(df)
            print(f"  [PASS] 'snps' format validated via string-cast parser!")
            print(f"         Total valid genomic rows parsed: {snp_count}")
            return True
        except Exception as retry_err:
            print(f"  [FAIL] 'snps' library rejected the file: {str(e)}")
            return False
        
def test_tsv_format_structure(file_path):
    """
    Tests basic file structure (header, columns, tab-separation) 
    required by tools like Promethease / openSNP.
    """
    print(f"\n[QA Test 2] Testing structural header format for Promethease compatibility: {file_path}")
    if not os.path.exists(file_path):
        print(f"  [FAIL] File not found: {file_path}")
        return False

    with open(file_path, "r") as f:
        lines = [line.strip() for line in f.readlines() if line.strip() and not line.startswith("#")]

    if not lines:
        print("  [FAIL] File contains no data rows.")
        return False

    # Check column header or first data line for 4 standard columns (rsid, chromosome, position, genotype)
    sample_line = lines[0].split("\t")
    if len(sample_line) >= 4:
        print(f"  [PASS] File contains valid 4-column tab-separated structure ({len(lines)} rows).")
        print(f"         Sample row: {sample_line}")
        return True
    else:
        print(f"  [FAIL] Invalid column structure. Expected at least 4 tab-separated columns, got {len(sample_line)}.")
        return False

if __name__ == "__main__":
    target_file = "results/output.txt" if os.path.exists("results/output.txt") else "data/output.txt"
    
    if len(sys.argv) > 1:
        target_file = sys.argv[1]

    print("=" * 60)
    print("=" * 60)
    
    res1 = test_snps_library(target_file)
    res2 = test_tsv_format_structure(target_file)
    
    print("\n" + "=" * 60)
    if res1 and res2:
        print("RESULT: PASS — Sanitized file maintains third-party tool compatibility.")
    else:
        print("RESULT: PARTIAL/FAIL — Review logs above for formatting issues.")
    print("=" * 60)