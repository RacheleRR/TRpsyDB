"""
Meta router
GET /api/meta/         → database stats
GET /api/meta/gene_sets → list all gene sets in db
"""
from fastapi import APIRouter
from backend.db.database import get_connection, get_db_stats

router = APIRouter()

@router.get("/")
def get_meta():
    return get_db_stats()

@router.get("/gene_sets")
def list_gene_sets():
    with get_connection() as conn:
        rows = conn.execute(
            "SELECT name, description, phenotype, n_genes, source_url FROM gene_sets"
        ).fetchall()
    return [dict(r) for r in rows]
