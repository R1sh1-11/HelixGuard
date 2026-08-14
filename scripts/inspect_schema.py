import sqlite3

conn = sqlite3.connect("data/blocklist.db")
print("blocklist sample:")
for r in conn.execute("SELECT rsid, GeneSymbol, ReferenceAllele, ReferenceAlleleVCF, AlternateAlleleVCF FROM blocklist LIMIT 5"):
    print("  ", r)

print("\nrsid types:", conn.execute("SELECT DISTINCT typeof(rsid) FROM blocklist").fetchall())
print("has rs prefix:", conn.execute("SELECT COUNT(*) FROM blocklist WHERE CAST(rsid AS TEXT) LIKE 'rs%'").fetchone()[0], "of 8942")

print("\nld_neighbors sample:")
for r in conn.execute("SELECT * FROM ld_neighbors LIMIT 5"):
    print("  ", r)
print("distinct targets:", conn.execute("SELECT COUNT(DISTINCT target_rsid) FROM ld_neighbors").fetchone()[0])
conn.close()