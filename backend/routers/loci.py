"""
Loci router — new schema (trs, tr_established, tr_regulatory, tr_brain_qtl)
GET /api/loci/              → paginated browse, filter by type/region/regulatory
GET /api/loci/search?q=HTT  → search by gene name or tr_id
GET /api/loci/{tr_id}       → full locus detail card
"""
from fastapi import APIRouter, HTTPException, Query
from typing import Optional
from backend.db.database import query, query_one

router = APIRouter()

# ── Browse ────────────────────────────────────────────────────────────────────
@router.get("/")
def list_loci(
    limit:       int            = Query(50, le=500),
    offset:      int            = Query(0),
    locus_type:  Optional[str]  = Query(None, description="mendelian|functional|clinvar|polymorphic"),
    gene_region: Optional[str]  = Query(None),
    in_promoter: Optional[int]  = Query(None, description="1 to filter"),
    in_brain_se: Optional[int]  = Query(None),
    in_ccre:     Optional[int]  = Query(None),
    has_qtl:     Optional[int]  = Query(None),
    rexprt_min:  Optional[float]= Query(None),
    motif_size:  Optional[str]  = Query(None, description="str|vntr"),
    sort:        Optional[str]  = Query("rexprt_desc"),
):
    clauses, params = [], []

    # locus type filter
    if locus_type == "mendelian":
        clauses.append("t.tr_id IN (SELECT DISTINCT tr_id FROM tr_established)")
    elif locus_type == "functional":
        clauses.append("t.tr_id IN (SELECT DISTINCT tr_id FROM tr_clinvar_functional WHERE locus_type='functional_VNTR')")
    elif locus_type == "clinvar":
        clauses.append("t.tr_id IN (SELECT DISTINCT tr_id FROM tr_clinvar_functional WHERE locus_type='ClinVar')")
    elif locus_type == "polymorphic":
        clauses.append("t.tr_id NOT IN (SELECT DISTINCT tr_id FROM tr_established)")
        clauses.append("t.tr_id NOT IN (SELECT DISTINCT tr_id FROM tr_clinvar_functional)")

    if gene_region:
        clauses.append("t.gene_region LIKE ?"); params.append(f"%{gene_region}%")
    if in_promoter == 1:
        clauses.append("r.in_promoter = 1")
    if in_brain_se == 1:
        clauses.append("r.in_brain_super_enhancer = 1")
    if in_ccre == 1:
        clauses.append("r.in_any_ccre = 1")
    if has_qtl == 1:
        clauses.append("t.tr_id IN (SELECT DISTINCT tr_id FROM tr_brain_qtl)")
    if rexprt_min is not None:
        clauses.append("rx.ensembleScore >= ?"); params.append(rexprt_min)
    if motif_size == "str":
        clauses.append("t.motif_size_bp <= 6")
    elif motif_size == "vntr":
        clauses.append("t.motif_size_bp >= 7")

    where = ("WHERE " + " AND ".join(clauses)) if clauses else ""

    sort_map = {
        "rexprt_desc": "rx.ensembleScore DESC NULLS LAST",
        "rexprt_asc":  "rx.ensembleScore ASC NULLS LAST",
        "gene_asc":    "t.gene_name ASC",
        "qtl_desc":    "(SELECT COUNT(*) FROM tr_brain_qtl q WHERE q.tr_id=t.tr_id) DESC",
    }
    order = sort_map.get(sort, "rx.ensembleScore DESC NULLS LAST")

    count_sql = f"""
        SELECT COUNT(*) as n FROM trs t
        LEFT JOIN tr_regulatory r  ON t.tr_id = r.tr_id
        LEFT JOIN tr_rexprt    rx  ON t.tr_id = rx.tr_id
        {where}
    """
    total = (query_one(count_sql, tuple(params)) or {}).get("n", 0)

    rows = query(f"""
        SELECT
            t.tr_id, t.chr, t.start, t.end,
            t.gene_name, t.gene_region, t.canonical_motif, t.motif_size_bp,
            t.n_sources, t.has_eVNTR,
            rx.ensembleScore as rexprt_score, rx.pLi, rx.loeuf,
            r.in_promoter, r.in_brain_super_enhancer, r.in_any_ccre,
            (SELECT COUNT(*) FROM tr_brain_qtl q WHERE q.tr_id=t.tr_id) as n_qtl,
            CASE WHEN e.tr_id IS NOT NULL THEN 'mendelian'
                 WHEN c.tr_id IS NOT NULL THEN c.locus_type
                 ELSE 'polymorphic' END as locus_type,
            e.disease_name
        FROM trs t
        LEFT JOIN tr_regulatory       r  ON t.tr_id = r.tr_id
        LEFT JOIN tr_rexprt           rx ON t.tr_id = rx.tr_id
        LEFT JOIN tr_established      e  ON t.tr_id = e.tr_id
        LEFT JOIN tr_clinvar_functional c ON t.tr_id = c.tr_id
        {where}
        GROUP BY t.tr_id
        ORDER BY {order}
        LIMIT ? OFFSET ?
    """, tuple(params) + (limit, offset))

    return {"total": total, "limit": limit, "offset": offset, "results": rows}

# ── Search ────────────────────────────────────────────────────────────────────
@router.get("/search")
def search_loci(q: str = Query(..., min_length=2)):
    rows = query("""
        SELECT
            t.tr_id, t.chr, t.start, t.end,
            t.gene_name, t.gene_region, t.canonical_motif,
            rx.ensembleScore as rexprt_score,
            CASE WHEN e.tr_id IS NOT NULL THEN 'mendelian'
                 WHEN c.tr_id IS NOT NULL THEN c.locus_type
                 ELSE 'polymorphic' END as locus_type,
            e.disease_name
        FROM trs t
        LEFT JOIN tr_rexprt           rx ON t.tr_id = rx.tr_id
        LEFT JOIN tr_established      e  ON t.tr_id = e.tr_id
        LEFT JOIN tr_clinvar_functional c ON t.tr_id = c.tr_id
        WHERE t.gene_name LIKE ?
           OR t.tr_id = ?
           OR e.disease_name LIKE ?
        GROUP BY t.tr_id
        ORDER BY rx.ensembleScore DESC NULLS LAST
        LIMIT 50
    """, (f"%{q}%", q, f"%{q}%"))
    return {"query": q, "n": len(rows), "results": rows}

# ── Full locus detail ─────────────────────────────────────────────────────────
@router.get("/{tr_id}")
def get_locus(tr_id: str):
    # core TR
    locus = query_one("SELECT * FROM trs WHERE tr_id = ?", (tr_id,))
    if not locus:
        raise HTTPException(status_code=404, detail=f"Locus '{tr_id}' not found")

    gene = locus.get("gene_name")

    # disease info
    established = query(
        "SELECT * FROM tr_established WHERE tr_id = ? ORDER BY id", (tr_id,)
    )

    # functional/clinvar
    clinvar = query(
        "SELECT * FROM tr_clinvar_functional WHERE tr_id = ?", (tr_id,)
    )

    # functional scores
    rexprt = query_one("SELECT * FROM tr_rexprt WHERE tr_id = ?", (tr_id,))

    # regulatory
    regulatory = query_one("SELECT * FROM tr_regulatory WHERE tr_id = ?", (tr_id,))

    # brain QTLs — top hits per cohort/qtl_type
    qtls = query("""
        SELECT cohort, tissue, qtl_type, gene_id, beta, p_value, q_value
        FROM tr_brain_qtl
        WHERE tr_id = ?
        ORDER BY p_value ASC
        LIMIT 50
    """, (tr_id,))

    # psychiatric literature
    psychiatric = query(
        "SELECT * FROM tr_psychiatric WHERE tr_id = ? ORDER BY evidence_tier", (tr_id,)
    )

    # gene-level psychiatric evidence
    gene_evidence = query(
        "SELECT * FROM gene_psychiatric_assoc WHERE gene_name = ? ORDER BY p_value",
        (gene,)
    ) if gene else []

    # PheWAS
    phewas = query(
        "SELECT * FROM tr_phenotype_assoc WHERE tr_id = ? ORDER BY p_value", (tr_id,)
    )

    # sc-eTR expression
    sc_etr = query(
        "SELECT * FROM tr_expression WHERE tr_id = ? ORDER BY p_fdr", (tr_id,)
    )

    # source catalog count
    n_qtl = (query_one(
        "SELECT COUNT(*) as n FROM tr_brain_qtl WHERE tr_id = ?", (tr_id,)
    ) or {}).get("n", 0)

    return {
        "locus":        locus,
        "established":  established,
        "clinvar":      clinvar,
        "rexprt":       rexprt,
        "regulatory":   regulatory,
        "qtls":         qtls,
        "n_qtl":        n_qtl,
        "psychiatric":  psychiatric,
        "gene_evidence":gene_evidence,
        "phewas":       phewas,
        "sc_etr":       sc_etr,
    }