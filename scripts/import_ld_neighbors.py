import sqlite3
import pandas as pd

df = pd.read_csv("data/ld_neighbors.csv")
conn = sqlite3.connect("data/blocklist.db")
df.to_sql("ld_neighbors", conn, if_exists="replace", index=False)
conn.execute("CREATE INDEX IF NOT EXISTS idx_target_rsid ON ld_neighbors(target_rsid)")
conn.commit()
conn.close()
print(f"Imported {len(df)} LD pairs")