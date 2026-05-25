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

