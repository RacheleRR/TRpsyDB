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
    if _is_turso():
        return libsql.connect(database=TURSO_URL, auth_token=TURSO_AUTH_TOKEN)
    if not DB_PATH.exists():
        raise FileNotFoundError(f"DB not found at {DB_PATH}.")
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
    def n(sql, p=()):
        return (query_one(sql, p) or {}).get("n", 0)
    return {
        "schema_version":      (query_one("SELECT value FROM db_meta WHERE key='schema_version'") or {}).get("value","1.0"),
        "n_loci":              n("SELECT COUNT(*) as n FROM trs"),
        "n_established":       n("SELECT COUNT(DISTINCT tr_id) as n FROM tr_established"),
        "n_clinvar_functional":n("SELECT COUNT(DISTINCT tr_id) as n FROM tr_clinvar_functional"),
        "n_with_rexprt":       n("SELECT COUNT(*) as n FROM tr_rexprt"),
        "n_with_qtl":          n("SELECT COUNT(DISTINCT tr_id) as n FROM tr_brain_qtl"),
        "n_in_promoter":       n("SELECT COUNT(*) as n FROM tr_regulatory WHERE in_promoter=1"),
        "n_in_brain_se":       n("SELECT COUNT(*) as n FROM tr_regulatory WHERE in_brain_super_enhancer=1"),
        "n_in_ccre":           n("SELECT COUNT(*) as n FROM tr_regulatory WHERE in_any_ccre=1"),
        "using_turso":         _is_turso(),
        "db_path":             str(TURSO_URL) if _is_turso() else str(DB_PATH),
    }

if __name__ == "__main__":
    import json
    print(json.dumps(get_db_stats(), indent=2))