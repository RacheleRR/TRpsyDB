"""
TRpsyDB — database connection
Supports local SQLite (dev) and Turso cloud (production on Render)

Local dev:   python runs against scripts/data/trs.db automatically
Production:  set TURSO_URL + TURSO_AUTH_TOKEN in Render environment vars
"""
import os
import sqlite3
from pathlib import Path

# ── Try to import libsql for Turso support ────────────────────────────────────
try:
    import libsql_experimental as libsql
    LIBSQL_AVAILABLE = True
except ImportError:
    LIBSQL_AVAILABLE = False

# ── Config ────────────────────────────────────────────────────────────────────
TURSO_URL        = os.getenv("TURSO_URL")
TURSO_AUTH_TOKEN = os.getenv("TURSO_AUTH_TOKEN")

# Local DB path — looks for trs.db relative to this file or via env var
_here    = Path(__file__).resolve().parent
DB_PATH  = Path(os.getenv("DB_PATH",
           str(_here.parent / "scripts" / "data" / "trs.db")))

def _is_turso() -> bool:
    return bool(TURSO_URL and TURSO_AUTH_TOKEN and LIBSQL_AVAILABLE)

# ── Connection ────────────────────────────────────────────────────────────────
def get_connection():
    """
    Returns a database connection.
    Turso cloud in production, local SQLite in development.
    """
    if _is_turso():
        conn = libsql.connect(
            database=TURSO_URL,
            auth_token=TURSO_AUTH_TOKEN,
        )
    else:
        if not DB_PATH.exists():
            raise FileNotFoundError(
                f"Database not found at {DB_PATH}. "
                f"Run 01_create_schema.py + 03_ingest_trexplorer.py first, "
                f"or set TURSO_URL + TURSO_AUTH_TOKEN for Turso."
            )
        conn = sqlite3.connect(str(DB_PATH))
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        conn.execute("PRAGMA journal_mode = WAL")
    return conn

# ── Query helpers ─────────────────────────────────────────────────────────────
def query(sql: str, params: tuple = ()) -> list[dict]:
    """Execute a SELECT and return list of dicts."""
    conn = get_connection()
    try:
        cur = conn.execute(sql, params)
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, row)) for row in cur.fetchall()]
    finally:
        conn.close()

def query_one(sql: str, params: tuple = ()) -> dict | None:
    """Execute a SELECT and return first row as dict or None."""
    results = query(sql, params)
    return results[0] if results else None

def execute(sql: str, params: tuple = ()) -> None:
    """Execute a write statement (INSERT/UPDATE/DELETE)."""
    conn = get_connection()
    try:
        conn.execute(sql, params)
        conn.commit()
    finally:
        conn.close()

# ── Stats for /api/meta ───────────────────────────────────────────────────────
def get_db_stats() -> dict:
    """Quick stats for the /api/meta endpoint."""
    return {
        "schema_version":     (query_one("SELECT value FROM db_meta WHERE key='schema_version'") or {}).get("value", "unknown"),
        "n_loci":             (query_one("SELECT COUNT(*) as n FROM trs") or {}).get("n", 0),
        "n_established":      (query_one("SELECT COUNT(*) as n FROM tr_established") or {}).get("n", 0),
        "n_clinvar_functional":(query_one("SELECT COUNT(*) as n FROM tr_clinvar_functional") or {}).get("n", 0),
        "n_with_rexprt":      (query_one("SELECT COUNT(*) as n FROM tr_rexprt") or {}).get("n", 0),
        "n_with_qtl":         (query_one("SELECT COUNT(*) as n FROM tr_brain_qtl") or {}).get("n", 0),
        "n_with_regulatory":  (query_one("SELECT COUNT(*) as n FROM tr_regulatory WHERE in_promoter=1 OR in_any_ccre=1") or {}).get("n", 0),
        "using_turso":        _is_turso(),
        "db_path":            str(TURSO_URL) if _is_turso() else str(DB_PATH),
    }

if __name__ == "__main__":
    import json
    print(json.dumps(get_db_stats(), indent=2))