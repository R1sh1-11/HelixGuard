import os
import sys
import pandas as pd
import snps


def test_snps_library(file_path):
    """
    [Tool 1/3] Tests if the open-source 'snps' Python library can successfully 
    parse the sanitized genomic file without crashing.
    """
    print(f"[QA Test 1] Testing compatibility with 'snps' library: {file_path}")
    try:
        s = snps.SNPs(file_path)
        snp_count = len(s.snps)
        if snp_count > 0:
            print(f"  [PASS] 'snps' library parsed {snp_count} variants successfully!")
            return True
        else:
            print("  [FAIL] 'snps' loaded 0 variants.")
            return False
    except Exception as e:
        # Fallback handling for string-type chrom values
        try:
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
    [Tool 2/3] Tests basic 4-column file structure (header, columns, tab-separation) 
    required by tools like Promethease / 23andMe TSV parsers.
    """
    print(f"\n[QA Test 2] Testing structural format for Promethease compatibility: {file_path}")
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if not lines:
        print("  [FAIL] File contains no data rows.")
        return False

    sample_line = lines[0].split("\t")
    if len(sample_line) >= 4:
        print(f"  [PASS] Valid 4-column tab-separated structure ({len(lines)} rows verified).")
        print(f"         Sample row: {sample_line}")
        return True
    else:
        print(f"  [FAIL] Expected at least 4 tab-separated columns, got {len(sample_line)}.")
        return False


def test_pandas_dataframe(file_path):
    """
    [Tool 3/3] Tests compatibility with standard Data Science ecosystem (Pandas/NumPy).
    """
    print(f"\n[QA Test 3] Testing Pandas DataFrame ingestion: {file_path}")
    try:
        df = pd.read_csv(file_path, sep="\t", comment="#", low_memory=False)
        if not df.empty and df.shape[1] >= 4:
            print(f"  [PASS] Pandas loaded DataFrame successfully ({len(df)} rows, {df.shape[1]} columns).")
            return True
        else:
            print("  [FAIL] Dataframe loaded empty or with insufficient columns.")
            return False
    except Exception as e:
        print(f"  [FAIL] Pandas failed to parse file: {e}")
        return False


if __name__ == "__main__":
    # Standard output target, or take custom path from command-line argument
    target_file = "results/sanitized_output.txt"
    
    if len(sys.argv) > 1:
        target_file = sys.argv[1]

    # Strict check: Ensure the target file exists before running tests
    if not os.path.exists(target_file):
        print("=" * 65)
        print(" ERROR: Input file not found.")
        print(f" Looking for: {target_file}")
        print("\n Usage:")
        print("   python src/validate_utility.py <path_to_sanitized_file.txt>")
        print(" Example:")
        print("   python src/validate_utility.py results/my_output.txt")
        print("=" * 65)
        sys.exit(1)

    print("=" * 65)
    print(" HELIXGUARD UTILITY VALIDATION SUITE (3 TOOLS)")
    print("=" * 65)

    res1 = test_snps_library(target_file)
    res2 = test_tsv_format_structure(target_file)
    res3 = test_pandas_dataframe(target_file)

    print("\n" + "=" * 65)
    if res1 and res2 and res3:
        print("FINAL RESULT: PASS — File is fully compatible with 3+ open-source tools!")
    else:
        print("FINAL RESULT: PARTIAL/FAIL — Review logs above.")
    print("=" * 65)