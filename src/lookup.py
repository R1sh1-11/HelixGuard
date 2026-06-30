import sqlite3

def is_flagged(rsid):
    rsid_int = int(rsid.replace("rs", ""))
    conn = sqlite3.connect("data/blocklist.db")
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM blocklist WHERE "RS# (dbSNP)" = ?', (rsid_int,))
    result = cursor.fetchone()
    conn.close()
    return result is not None, result

if __name__ == "__main__":
    flagged, data = is_flagged("rs429358")
    print(flagged, data)