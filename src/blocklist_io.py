"""
Central blocklist access. The column was renamed from "RS# (dbSNP)" to "rsid"
by Kritha's week 9 rebuild. SQLite silently treats an unknown double-quoted
identifier as a string literal instead of erroring, so stale column names
return garbage. All blocklist reads go through here.
"""

import sqlite3

DB_PATH = "data/blocklist.db"


def _rsid_column(conn):
    cols = [d[1] for d in conn.execute("PRAGMA table_info(blocklist)")]
    for candidate in ("rsid", "RS# (dbSNP)"):
        if candidate in cols:
            return candidate
    raise RuntimeError(f"No rsID column found. Columns: {cols}")


def load_reference_map(db_path=DB_PATH):
    """rsID string -> homozygous reference genotype string."""
    conn = sqlite3.connect(db_path)
    col = _rsid_column(conn)
    rows = conn.execute(
        f'SELECT "{col}", ReferenceAllele FROM blocklist'
    ).fetchall()
    conn.close()
    return {
        str(r[0]).strip(): str(r[1]).strip().upper()
        for r in rows if r[0] is not None and r[1]
    }


def load_blocklist_rsids(db_path=DB_PATH):
    """Set of blocklisted rsID strings."""
    return set(load_reference_map(db_path).keys())