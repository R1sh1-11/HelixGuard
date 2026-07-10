import sqlite3

def get_ld_neighbors(rsid, r2_threshold=0.8):
    """
    Given a flagged rsID string (e.g. 'rs429358'), queries the precomputed
    SQLite index to instantly find correlated neighborhood SNPs.
    """
    # Standardize input formatting (e.g., '  RS429358 ' -> 'rs429358')
    rsid_clean = str(rsid).strip().lower()
    if not rsid_clean.startswith("rs"):
        rsid_clean = f"rs{rsid_clean}"

    conn = sqlite3.connect("data/blocklist.db")
    cursor = conn.cursor()
    
    # Check if our table exists yet
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='ld_neighbors'")
    if not cursor.fetchone():
        conn.close()
        return []
    
    # Query our fast pre-indexed table
    query = "SELECT neighbor_rsid FROM ld_neighbors WHERE target_rsid = ? AND r2 >= ?"
    cursor.execute(query, (rsid_clean, r2_threshold))
    rows = cursor.fetchall()
    conn.close()
    
    # Return a clean list of neighbor rsID strings
    return [row[0] for row in rows]

if __name__ == "__main__":
    # Quick local test verification
    example_rsid = "rs429358" 
    results = get_ld_neighbors(example_rsid)
    print(f"\n[Test] Found {len(results)} highly linked neighbors for {example_rsid}:")
    print(results)