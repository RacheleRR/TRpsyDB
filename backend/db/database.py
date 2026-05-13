"""
psyTRdb — database connection and initialization
"""
import sqlite3
import os
from pathlib import Path

# Path to the SQLite database file
DB_DIR  = Path(__file__).parent
DB_PATH = DB_DIR / "psytrdb.sqlite"
SCHEMA_PATH = DB_DIR / "schema.sql"


def get_connection() -> sqlite3.Connection:
    """Return a database connection with row_factory so results come back as dicts."""
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row      # access columns by name: row["gene_symbol"]
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA journal_mode = WAL")   # better concurrent read performance
    return conn


def init_db() -> None:
    """Create all tables if they don't exist yet. Safe to run on every startup."""
    with get_connection() as conn:
        sql = SCHEMA_PATH.read_text()
        conn.executescript(sql)
    print(f"[DB] Initialized at {DB_PATH}")


def get_db_stats() -> dict:
    """Quick stats for the /api/meta endpoint."""
    with get_connection() as conn:
        n_loci    = conn.execute("SELECT COUNT(*) FROM loci").fetchone()[0]
        n_studies = conn.execute("SELECT COUNT(*) FROM studies").fetchone()[0]
        n_tier1   = conn.execute("SELECT COUNT(*) FROM loci WHERE evidence_tier = 1").fetchone()[0]
        n_tier2   = conn.execute("SELECT COUNT(*) FROM loci WHERE evidence_tier = 2").fetchone()[0]
        version   = conn.execute("SELECT value FROM db_meta WHERE key='schema_version'").fetchone()[0]
    return {
        "schema_version": version,
        "n_loci": n_loci,
        "n_studies": n_studies,
        "n_tier1_established": n_tier1,
        "n_tier2_strong": n_tier2,
    }


if __name__ == "__main__":
    init_db()
    print(get_db_stats())
