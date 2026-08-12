import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sqlite3
from src.blocklist_io import load_blocklist_rsids

DB = "data/blocklist.db"

if not os.path.exists(DB):
    print("ERROR: data/blocklist.db not found")
    raise SystemExit(1)

conn = sqlite3.connect(DB)

tables = [r[0] for r in conn.execute(
    "SELECT name FROM sqlite_master WHERE type='table'")]
print("tables:", tables)

bl_rows = conn.execute("SELECT COUNT(*) FROM blocklist").fetchone()[0]
cols = [d[1] for d in conn.execute("PRAGMA table_info(blocklist)")]
print(f"\nblocklist rows: {bl_rows:,}")
print("blocklist columns:", cols)
print("ReferenceAlleleVCF column present:", "ReferenceAlleleVCF" in cols)

print("\nper-gene reference allele spread:")
gene_variety = {}
for gene, allele, n in conn.execute(
        "SELECT GeneSymbol, ReferenceAllele, COUNT(*) FROM blocklist "
        "GROUP BY GeneSymbol, ReferenceAllele ORDER BY GeneSymbol"):
    print(f"  {gene:8} {allele:4} {n:,}")
    gene_variety[gene] = gene_variety.get(gene, 0) + 1

max_variety = max(gene_variety.values()) if gene_variety else 0
ref_fixed = max_variety > 1

if "ld_neighbors" in tables:
    ld_rows = conn.execute("SELECT COUNT(*) FROM ld_neighbors").fetchone()[0]
    ld_targets = {r[0] for r in conn.execute(
        "SELECT DISTINCT target_rsid FROM ld_neighbors")}
else:
    ld_rows = 0
    ld_targets = set()
print(f"\nld_neighbors pairs: {ld_rows:,}")

bl_rsids = load_blocklist_rsids()
overlap = len(ld_targets & bl_rsids)
pct = (overlap / len(ld_targets) * 100) if ld_targets else 0
print(f"ld targets still in blocklist: {overlap} of {len(ld_targets)} ({pct:.1f}%)")

conn.close()

genomes = [f for f in os.listdir("data")
           if f.startswith("genome_") and f.endswith(".txt")]
print(f"\ngenome files in data/: {len(genomes)}")
for g in sorted(genomes):
    print("  ", g)

print("\n" + "=" * 60)
print("VERDICTS")
print("=" * 60)
print(f"1. Reference alleles fixed : {'YES' if ref_fixed else 'NO -- still one per gene'}")
print(f"2. LD table populated      : {'YES' if ld_rows > 100 else f'NO -- only {ld_rows} pairs'}")
print(f"3. LD data still valid     : {'YES' if pct >= 50 else f'NO -- only {pct:.0f}% overlap, refetch needed'}")
print(f"4. Genome count            : {len(genomes)} of 10 {'OK' if len(genomes) >= 10 else '-- MISSING FILES'}")
print("=" * 60)