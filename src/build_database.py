import sqlite3
import pandas as pd

# Load blocklist into SQLite
df = pd.read_csv("data/blocklist.csv")

conn = sqlite3.connect("data/blocklist.db")
df.to_sql("blocklist", conn, if_exists="replace", index=False)
conn.close()
print("Database created.")

# Test the lookup
conn = sqlite3.connect("data/blocklist.db")
cursor = conn.cursor()
cursor.execute("SELECT * FROM blocklist WHERE \"RS# (dbSNP)\" = 429358")
print(cursor.fetchall())
conn.close()