import sqlite3
import pandas as pd

conn = sqlite3.connect("data/blocklist.db")
df = pd.read_sql("SELECT target_rsid, neighbor_rsid, r2 FROM ld_neighbors", conn)
conn.close()
df.to_csv("data/ld_neighbors.csv", index=False)
print(f"Exported {len(df)} LD pairs to data/ld_neighbors.csv")