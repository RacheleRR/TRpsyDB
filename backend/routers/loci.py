"""
Loci router — Mode 1 endpoints
GET /api/loci/search?q=RFC1              → search by gene name or locus_id
GET /api/loci/{locus_id}                 → full locus report card
GET /api/loci/{locus_id}/evidence        → all studies that tested this locus
GET /api/loci/gene/{gene_symbol}         → all TR loci within/near a gene
GET /api/loci/                           → browse all loci (paginated, filterable by tier)
"""
from fastapi import APIRouter, HTTPException, Query
from typing import Optional
from backend.db.database import get_connection

router = APIRouter()


def row_to_dict(row) -> dict:
    """Convert sqlite3.Row to plain dict."""
    return dict(row) if row else None


# ── Browse all loci ──────────────────────────────────────────
@router.get("/")
def list_loci(
    tier: Optional[int]  = Query(None, description="Filter by evidence tier (1-5)"),
    phenotype: Optional[str] = Query(None, description="Filter by phenotype e.g. SCZ"),
    limit: int = Query(50, le=500),
    offset: int = Query(0),
):
    with get_connection() as conn:
        clauses, params = [], []

        if tier is not None:
            clauses.append("evidence_tier = ?")
            params.append(tier)
        if phenotype:
            clauses.append("phenotypes LIKE ?")
            params.append(f"%{phenotype}%")

        where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
        total = conn.execute(
            f"SELECT COUNT(*) FROM loci {where}", params
        ).fetchone()[0]

        rows = conn.execute(
            f"""SELECT locus_id, gene_symbol, chrom, start, end, motif,
                       evidence_tier, phenotypes, rexprt_score, common_name
                FROM loci {where}
                ORDER BY evidence_tier, gene_symbol
                LIMIT ? OFFSET ?""",
            params + [limit, offset],
        ).fetchall()

    return {
        "total": total,
        "limit": limit,
        "offset": offset,
        "results": [row_to_dict(r) for r in rows],
    }


# ── Search by gene name or locus_id ─────────────────────────
@router.get("/search")
def search_loci(q: str = Query(..., min_length=2)):
    """
    Search by gene symbol (partial match) or exact locus_id.
    """
    with get_connection() as conn:
        rows = conn.execute(
            """SELECT locus_id, gene_symbol, chrom, start, end, motif,
                      evidence_tier, phenotypes, rexprt_score, common_name
               FROM loci
               WHERE gene_symbol LIKE ?
                  OR locus_id = ?
                  OR common_name LIKE ?
               ORDER BY evidence_tier, gene_symbol
               LIMIT 50""",
            (f"%{q}%", q, f"%{q}%"),
        ).fetchall()

    return {"query": q, "n": len(rows), "results": [row_to_dict(r) for r in rows]}


# ── Full locus report card ───────────────────────────────────
@router.get("/{locus_id:path}")
def get_locus(locus_id: str):
    """
    Full report card for a single locus.
    locus_id format: chrN:start-end   e.g. chr3:63912684-63912714
    """
    with get_connection() as conn:
        locus = conn.execute(
            "SELECT * FROM loci WHERE locus_id = ?", (locus_id,)
        ).fetchone()

        if not locus:
            raise HTTPException(status_code=404, detail=f"Locus '{locus_id}' not found")

        evidence = conn.execute(
            """SELECT le.*, s.title, s.first_author, s.year, s.journal,
                      s.pmid, s.cohort_size_cases, s.cohort_size_controls,
                      s.phenotype, s.sequencing_type
               FROM locus_evidence le
               JOIN studies s ON le.study_id = s.id
               WHERE le.locus_id = ?
               ORDER BY s.year DESC""",
            (locus_id,),
        ).fetchall()

    return {
        "locus": row_to_dict(locus),
        "evidence": [row_to_dict(e) for e in evidence],
    }


# ── Full locus report card ───────────────────────────────────
@router.get("/{locus_id:path}")
def get_locus(locus_id: str):
    """
    Full report card for a single locus.
    locus_id format: chrN:start-end   e.g. chr3:63912684-63912714
    """
    with get_connection() as conn:
        locus = conn.execute(
            "SELECT * FROM loci WHERE locus_id = ?", (locus_id,)
        ).fetchone()

        if not locus:
            raise HTTPException(status_code=404, detail=f"Locus '{locus_id}' not found")

        evidence = conn.execute(
            """SELECT le.*, s.title, s.first_author, s.year, s.journal,
                      s.pmid, s.cohort_size_cases, s.cohort_size_controls,
                      s.phenotype, s.sequencing_type
               FROM locus_evidence le
               JOIN studies s ON le.study_id = s.id
               WHERE le.locus_id = ?
               ORDER BY s.year DESC""",
            (locus_id,),
        ).fetchall()

    return {
        "locus": row_to_dict(locus),
        "evidence": [row_to_dict(e) for e in evidence],
    }


# ── All loci near a gene ─────────────────────────────────────
@router.get("/gene/{gene_symbol}")
def get_loci_by_gene(gene_symbol: str):
    """All TR loci associated with a given gene symbol (case-insensitive)."""
    with get_connection() as conn:
        rows = conn.execute(
            """SELECT * FROM loci
               WHERE UPPER(gene_symbol) = UPPER(?)
               ORDER BY evidence_tier, tss_distance""",
            (gene_symbol,),
        ).fetchall()

    if not rows:
        raise HTTPException(
            status_code=404,
            detail=f"No TR loci found for gene '{gene_symbol}'. "
                   "Check spelling or try a search.",
        )

    return {"gene": gene_symbol.upper(), "n": len(rows), "loci": [row_to_dict(r) for r in rows]}


