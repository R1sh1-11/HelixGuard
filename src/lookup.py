import sqlite3

def is_flagged(rsid):
    try:
        #string cleaning logic to strip 'rs' and handle integers
        rsid_int = int(rsid.lower().replace("rs", "").strip())
    except ValueError:
        return False, None
        
    conn = sqlite3.connect("data/blocklist.db")
    cursor = conn.cursor()
    
    #select our new ReferenceAllele column instead of the entire row
    cursor.execute('SELECT "ReferenceAllele" FROM blocklist WHERE "RS# (dbSNP)" = ?', (rsid_int,))
    result = cursor.fetchone()
    conn.close()
    
    if result:
        #result is returned as a single-element tuple like ('TT',), extract string at index 0
        return True, result[0]
        
    return False, None

if __name__ == "__main__":
    #internal test to ensure it works end-to-end
    flagged, safe_genotype = is_flagged("rs429358")
    print(f"Is Flagged: {flagged} | Safe Replacement Allele: {safe_genotype}")