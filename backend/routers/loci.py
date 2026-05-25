"""
Loci router — new schema
GET /api/loci/              → paginated browse with filters
GET /api/loci/search?q=HTT  → search by gene or tr_id
GET /api/loci/{tr_id}       → full locus detail (batched + cached)
"""
from fastapi import APIRouter, HTTPException, Query
from typing import Optional
from functools import lru_cache
import json
from backend.db.database import query, query_one

router = APIRouter()

# ── FIX 2: cache locus detail results ────────────────────────────────────────
@lru_cache(maxsize=500)
def _get_locus_cached(tr_id: str) -> str:
    """Cache full locus detail as JSON string. 500 loci in RAM."""

    # FIX 1: batch core data into ONE query instead of 8 round trips
    locus = query_one("""
        SELECT
            t.*,
            rx.ensembleScore, rx.ensembleBinary, rx.ensembleMax,
            rx.SVM, rx.XGB, rx.pLi, rx.loeuf,
            rx.gc_content, rx.per_g, rx.per_c, rx.per_a, rx.per_t,
            rx.eSTR, rx.opReg, rx.promoter, rx.HG19_ID,
            r.in_promoter, r.promoter_id,
            r.ccre_class, r.ccre_id, r.in_any_ccre, r.in_brain_ccre,
            r.in_any_super_enhancer, r.in_brain_super_enhancer,
            r.n_super_enhancer_tissues, r.super_enhancer_tissues,
            r.SE_gene_overlap, r.SE_gene_closest, r.SE_gene_ABC, r.SE_cas_value,
            r.in_any_enhancer, r.in_brain_enhancer, r.n_enhancer_tissues,
            r.dist_tss_bp, r.nearest_tss_gene, r.dist_splice_bp
        FROM trs t
        LEFT JOIN tr_rexprt     rx ON t.tr_id = rx.tr_id
        LEFT JOIN tr_regulatory r  ON t.tr_id = r.tr_id
        WHERE t.tr_id = ?
    """, (tr_id,))

    if not locus:
        return json.dumps(None)

    gene = locus.get("gene_name")

    # disease + psychiatric + gene evidence in one batch
    established  = query("SELECT * FROM tr_established WHERE tr_id=? ORDER BY id", (tr_id,))
    clinvar      = query("SELECT * FROM tr_clinvar_functional WHERE tr_id=?", (tr_id,))
    psychiatric  = query("SELECT * FROM tr_psychiatric WHERE tr_id=? ORDER BY evidence_tier", (tr_id,))
    gene_evid    = query("SELECT * FROM gene_psychiatric_assoc WHERE gene_name=? ORDER BY p_value", (gene,)) if gene else []
    phewas       = query("SELECT * FROM tr_phenotype_assoc WHERE tr_id=? ORDER BY p_value", (tr_id,))
    sc_etr       = query("SELECT * FROM tr_expression WHERE tr_id=? ORDER BY p_fdr", (tr_id,))

    # QTL summary — top 20 only for initial load (rest loaded lazily)
    qtls = query("""
        SELECT cohort, tissue, qtl_type, gene_id, beta, p_value, q_value
        FROM tr_brain_qtl WHERE tr_id=?
        ORDER BY p_value ASC LIMIT 20
    """, (tr_id,))

    n_qtl = (query_one("SELECT COUNT(*) as n FROM tr_brain_qtl WHERE tr_id=?", (tr_id,)) or {}).get("n", 0)

    result = {
        "locus":        locus,
        "established":  established,
        "clinvar":      clinvar,
        "psychiatric":  psychiatric,
        "gene_evidence":gene_evid,
        "phewas":       phewas,
        "sc_etr":       sc_etr,
        "qtls":         qtls,
        "n_qtl":        n_qtl,
    }
    return json.dumps(result, default=str)

# ── Browse ────────────────────────────────────────────────────────────────────
@router.get("/")
def list_loci(
    limit:       int            = Query(50,   le=500,  description="Max results per page"),
    offset:      int            = Query(0,            description="Pagination offset"),
    locus_type:  Optional[str]  = Query(None,         description="mendelian|functional|clinvar|polymorphic"),
    gene_region: Optional[str]  = Query(None,         description="Partial match e.g. Coding, Intronic, UTR, Promoter, Intergenic"),
    in_promoter: Optional[int]  = Query(None,         description="1 = only TRs overlapping an EPD promoter"),
    in_brain_se: Optional[int]  = Query(None,         description="1 = only TRs in a brain super enhancer (SEdb)"),
    in_ccre:     Optional[int]  = Query(None,         description="1 = only TRs overlapping an ENCODE cCRE"),
    has_qtl:     Optional[int]  = Query(None,         description="1 = only TRs with at least one brain QTL hit"),
    has_evntr:   Optional[int]  = Query(None,         description="1 = only TRs flagged as eVNTR (Garg 2021)"),
    rexprt_min:  Optional[float]= Query(None, ge=0, le=1, description="Minimum RExPRT ensemble score (0–1)"),
    motif_size:  Optional[str]  = Query(None,         description="str = motif ≤6 bp | vntr = motif ≥7 bp"),
    sort:        Optional[str]  = Query("rexprt_desc",description="rexprt_desc|rexprt_asc|gene_asc|qtl_desc|sources_desc"),
    q:           Optional[str]  = Query(None,         description="Free-text search: gene name, tr_id, or disease name"),
):
    clauses, params = [], []

    if q:
        clauses.append("(t.gene_name LIKE ? OR t.tr_id = ? OR e.disease_name LIKE ?)")
        params += [f"%{q}%", q, f"%{q}%"]

    if locus_type == "mendelian":
        clauses.append("e.tr_id IS NOT NULL")
    elif locus_type == "functional":
        clauses.append("c.locus_type = 'functional_VNTR'")
    elif locus_type == "clinvar":
        clauses.append("c.locus_type = 'ClinVar'")
    elif locus_type == "polymorphic":
        clauses.append("e.tr_id IS NULL AND c.tr_id IS NULL")

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
    if has_evntr == 1:
        clauses.append("t.has_eVNTR = 1")
    if rexprt_min is not None:
        clauses.append("rx.ensembleScore >= ?"); params.append(rexprt_min)
    if motif_size == "str":
        clauses.append("t.motif_size_bp <= 6")
    elif motif_size == "vntr":
        clauses.append("t.motif_size_bp >= 7")

    where = ("WHERE " + " AND ".join(clauses)) if clauses else ""

    sort_map = {
        "rexprt_desc": "rx.ensembleScore DESC",
        "rexprt_asc":  "rx.ensembleScore ASC",
        "gene_asc":    "t.gene_name ASC",
        "qtl_desc":    "n_qtl DESC",
        "sources_desc":"t.n_sources DESC",
    }
    order = sort_map.get(sort, "rx.ensembleScore DESC")

    joins = """
        FROM trs t
        LEFT JOIN tr_rexprt           rx ON t.tr_id = rx.tr_id
        LEFT JOIN tr_regulatory        r ON t.tr_id = r.tr_id
        LEFT JOIN tr_established       e ON t.tr_id = e.tr_id
        LEFT JOIN tr_clinvar_functional c ON t.tr_id = c.tr_id
    """

    total = (query_one(
        f"SELECT COUNT(DISTINCT t.tr_id) as n {joins} {where}", tuple(params)
    ) or {}).get("n", 0)

    rows = query(f"""
        SELECT
            t.tr_id, t.chr, t.start, t.end,
            t.gene_name, t.gene_region, t.canonical_motif, t.motif_size_bp,
            t.n_sources, t.has_eVNTR,
            rx.ensembleScore as rexprt_score,
            r.in_promoter, r.in_brain_super_enhancer, r.in_any_ccre,
            CASE WHEN e.tr_id IS NOT NULL THEN 'mendelian'
                 WHEN c.locus_type = 'functional_VNTR' THEN 'functional'
                 WHEN c.locus_type = 'ClinVar' THEN 'clinvar'
                 ELSE 'polymorphic' END as locus_type,
            e.disease_name,
            (SELECT COUNT(*) FROM tr_brain_qtl q WHERE q.tr_id=t.tr_id) as n_qtl
        {joins}
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
            e.disease_name,
            (SELECT COUNT(*) FROM tr_brain_qtl bq WHERE bq.tr_id=t.tr_id) as n_qtl
        FROM trs t
        LEFT JOIN tr_rexprt            rx ON t.tr_id = rx.tr_id
        LEFT JOIN tr_established        e ON t.tr_id = e.tr_id
        LEFT JOIN tr_clinvar_functional c ON t.tr_id = c.tr_id
        WHERE t.gene_name LIKE ?
           OR t.tr_id = ?
           OR e.disease_name LIKE ?
           OR t.canonical_motif = ?
        GROUP BY t.tr_id
        ORDER BY
            CASE WHEN e.tr_id IS NOT NULL THEN 0 ELSE 1 END,
            rx.ensembleScore DESC
        LIMIT 50
    """, (f"%{q}%", q, f"%{q}%", q.upper()))
    return {"query": q, "n": len(rows), "results": rows}

# ── Full locus detail (cached + batched) ─────────────────────────────────────
@router.get("/{tr_id}")
def get_locus(tr_id: str):
    cached = _get_locus_cached(tr_id)
    result = json.loads(cached)
    if result is None:
        raise HTTPException(status_code=404, detail=f"Locus '{tr_id}' not found")
    return result

# ── QTLs lazy load endpoint (FIX 4) ──────────────────────────────────────────
@router.get("/{tr_id}/qtls")
def get_qtls(tr_id: str, limit: int = Query(100, le=500)):
    """Lazy-loaded QTL endpoint — called after page renders."""
    rows = query("""
        SELECT cohort, tissue, qtl_type, gene_id, beta, se, p_value, q_value
        FROM tr_brain_qtl WHERE tr_id=?
        ORDER BY p_value ASC
        LIMIT ?
    """, (tr_id, limit))
    return {"tr_id": tr_id, "n": len(rows), "qtls": rows}