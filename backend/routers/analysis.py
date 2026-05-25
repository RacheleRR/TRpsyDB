"""
Analysis router — Mode 2 endpoints
POST /api/analysis/enrich     → submit a TR list, get back enrichment results
GET  /api/analysis/{job_id}   → check job status / retrieve results
"""
import uuid
import json
import math
import requests
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Optional
from scipy import stats

from backend.db.database import query, query_one, execute

router = APIRouter()


# ── Input schema ─────────────────────────────────────────────
class TRRecord(BaseModel):
    chrom: str                          # e.g. chr3
    start: int
    end: int
    motif: Optional[str] = None
    p_value: Optional[float] = None


class EnrichmentRequest(BaseModel):
    trs: List[TRRecord]
    label: Optional[str] = "My TR list"


# ── Helper: overlap with database loci ───────────────────────
#? should we do only matchy matchy or also overlappy ?
def find_database_matches(trs: List[TRRecord]) -> list:
    all_loci = query("""
        SELECT t.tr_id, t.chr, t.start, t.end, t.gene_name,
               rx.ensembleScore as rexprt_score,
               CASE WHEN e.tr_id IS NOT NULL THEN 1 ELSE 0 END as is_established,
               e.disease_name
        FROM trs t
        LEFT JOIN tr_rexprt      rx ON t.tr_id = rx.tr_id
        LEFT JOIN tr_established e  ON t.tr_id = e.tr_id
        GROUP BY t.tr_id
    """)

    matches = []
    for tr in trs:
        hit = None
        for locus in all_loci:
            if (locus["chr"] == tr.chrom and
                    tr.start < locus["end"] and
                    tr.end > locus["start"]):
                hit = dict(locus)
                break
        matches.append({
            "input":    tr.dict(),
            "db_match": hit,
            "novel":    hit is None,
        })
    return matches

# ── Helper: gene set enrichment (Fisher's exact) ─────────────
def gene_set_enrichment(input_genes: set) -> list:
    genome_gene_count = (query_one(
        "SELECT COUNT(DISTINCT gene_name) as n FROM trs WHERE gene_name IS NOT NULL"
    ) or {}).get("n", 20000)

    # build gene sets from gene_psychiatric_assoc
    sources = query(
        "SELECT DISTINCT source, phenotype FROM gene_psychiatric_assoc"
    )
    results = []
    for src in sources:
        gs_genes = set(r["gene_name"] for r in query(
            "SELECT gene_name FROM gene_psychiatric_assoc WHERE source=? AND phenotype=?",
            (src["source"], src["phenotype"])
        ))
        overlap = input_genes & gs_genes
        a = len(overlap)
        b = len(input_genes) - a
        c = len(gs_genes) - a
        d = genome_gene_count - a - b - c
        if a == 0:
            results.append({
                "gene_set":     f"{src['source']} ({src['phenotype']})",
                "phenotype":    src["phenotype"],
                "n_overlap":    0,
                "overlap_genes":[],
                "odds_ratio":   0,
                "p_value":      1.0,
                "significant":  False,
            })
            continue
        _, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        or_val = (a * d) / (b * c) if (b * c) > 0 else float("inf")
        results.append({
            "gene_set":     f"{src['source']} ({src['phenotype']})",
            "phenotype":    src["phenotype"],
            "n_overlap":    a,
            "overlap_genes":sorted(overlap),
            "odds_ratio":   round(or_val, 3),
            "p_value":      round(p_value, 6),
            "significant":  p_value < 0.05,
        })

    # Sort by p-value
    results.sort(key=lambda x: x["p_value"])
    return results


# ── Helper: TSS / splice proximity summary ───────────────────
# def proximity_summary(trs: List[TRRecord], conn) -> dict:
#     """
#     For each input TR, annotate TSS and splice junction proximity
#     by matching to nearest locus in database (or a pre-computed BED).
#     Returns a summary dict.
#     """
#     # For now: count how many input TRs are within 1kb of a TSS
#     # This will be extended once we have a full TSS BED file loaded
#     matched = find_database_matches(trs, conn)
#     known_with_tss = [
#         m["db_match"]["locus_id"]
#         for m in matched
#         if m["db_match"] and m["db_match"].get("tss_distance") is not None
#     ]
#     return {
#         "n_with_tss_annotation": len(known_with_tss),
#         "note": "Full TSS proximity for novel TRs requires genome annotation file (coming soon)"
#     }


# ── Helper: ORA via g:Profiler public API ────────────────────
GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP", "HP"]

SOURCE_LABELS = {
    "GO:BP": "GO Biological Process",
    "GO:MF": "GO Molecular Function",
    "GO:CC": "GO Cellular Component",
    "KEGG":  "KEGG",
    "REAC":  "Reactome",
    "WP":    "WikiPathways",
    "HP":    "Human Phenotype Ontology",
}

SOURCE_COLORS = {
    "GO:BP": "#2563a8",
    "GO:MF": "#1a6b3c",
    "GO:CC": "#7a4f1a",
    "KEGG":  "#8b3a8b",
    "REAC":  "#c0392b",
    "WP":    "#e67e22",
    "HP":    "#2c7a7a",
}

def run_gprofiler_ora(gene_list: list) -> dict:
    """
    Run ORA via g:Profiler public REST API.
    Returns structured results grouped by source, ready for frontend rendering.
    Gracefully handles API errors — analysis continues without ORA if it fails.
    """
    if not gene_list:
        return {"status": "skipped", "reason": "No genes to analyse", "results": []}

    try:
        response = requests.post(
            GPROFILER_URL,
            json={
                "organism":                     "hsapiens",
                "query":                        gene_list,
                "sources":                      SOURCES,
                "significance_threshold_method": "fdr",
                "user_threshold":               0.05,
                "no_evidences":                 True,   # faster, we don't need evidence codes
                "combined":                     False,
                "measure_underrepresentation":  False,
                "domain_scope":                 "annotated",
            },
            timeout=20,
        )
        response.raise_for_status()
        data = response.json()

    except requests.exceptions.Timeout:
        return {"status": "error", "reason": "g:Profiler API timed out", "results": []}
    except requests.exceptions.RequestException as e:
        return {"status": "error", "reason": f"g:Profiler API error: {str(e)}", "results": []}

    raw_results = data.get("result", [])

    if not raw_results:
        return {"status": "ok", "n_significant": 0, "results": []}

    # Parse and structure results
    parsed = []
    for r in raw_results:
        source = r.get("source", "")
        p_val  = r.get("p_value", 1.0)
        parsed.append({
            "term_id":          r.get("native", ""),
            "term_name":        r.get("name", ""),
            "source":           source,
            "source_label":     SOURCE_LABELS.get(source, source),
            "source_color":     SOURCE_COLORS.get(source, "#888"),
            "p_value":          round(p_val, 8),
            "p_value_log10":    round(-math.log10(p_val + 1e-300), 3),
            "fdr":              round(r.get("p_value", 1.0), 8),  # g:Profiler returns FDR-adjusted p
            "term_size":        r.get("term_size", 0),
            "query_size":       r.get("query_size", 0),
            "intersection_size": r.get("intersection_size", 0),
            "intersection_genes": r.get("intersections", []),
            "precision":        round(r.get("precision", 0), 4),
            "recall":           round(r.get("recall", 0), 4),
        })

    # Sort by p-value
    parsed.sort(key=lambda x: x["p_value"])

    # Group by source for frontend rendering
    by_source = {}
    for item in parsed:
        src = item["source"]
        if src not in by_source:
            by_source[src] = {
                "source":       src,
                "source_label": SOURCE_LABELS.get(src, src),
                "source_color": SOURCE_COLORS.get(src, "#888"),
                "terms":        [],
            }
        by_source[src]["terms"].append(item)

    return {
        "status":        "ok",
        "n_significant": len(parsed),
        "n_genes_query": len(gene_list),
        "results":       parsed[:100],      # cap at 100 for payload size
        "by_source":     list(by_source.values()),
        "top_terms":     parsed[:10],       # top 10 for summary display
    }


# ── Main endpoint ─────────────────────────────────────────────
@router.post("/enrich")
def submit_enrichment(req: EnrichmentRequest):
    """
    Submit a list of TR coordinates for enrichment analysis.
    Returns results synchronously for now (fast enough for <10k TRs).
    """
    if len(req.trs) > 10_000:
        raise HTTPException(status_code=400, detail="Maximum 10,000 TRs per submission.")
    if len(req.trs) == 0:
        raise HTTPException(status_code=400, detail="TR list is empty.")

    job_id = str(uuid.uuid4())
    # 1. Match input TRs to database loci
    matches = find_database_matches(req.trs)
    # 2. Gene set enrichment
    input_genes = set(
        m["db_match"]["gene_name"]
        for m in matches
        if m["db_match"] and m["db_match"].get("gene_name")
    )

    enrichment = gene_set_enrichment(input_genes) if input_genes else []

    # # 3. Proximity summary
    # proximity = proximity_summary(req.trs, conn)
    # 4. ORA via g:Profiler
    ora        = run_gprofiler_ora(sorted(input_genes))
    # 5. Novel vs known summary
    n_novel   = sum(1 for m in matches if m["novel"])
    n_known   = len(matches) - n_novel
    n_estab   = sum(1 for m in matches if m["db_match"] and m["db_match"].get("is_established"))


    result = {
        "job_id": job_id, "label": req.label, "n_input": len(req.trs),
        "summary": {
            "n_known_in_db":    n_known,
            "n_novel":          n_novel,
            "n_established":    n_estab,
            "n_unique_genes":   len(input_genes),
        },
        "matches":          matches,
        "gene_set_enrichment": enrichment,
        "ora":              ora,
    }

    try:
        execute("""
            INSERT INTO analysis_jobs (id, status, input_n_trs, completed_at, result_json)
            VALUES (?, 'done', ?, datetime('now'), ?)
        """, (job_id, len(req.trs), json.dumps(result)))
    except Exception:
        pass  # analysis_jobs table may not exist — non-fatal

    return result


@router.get("/{job_id}")
def get_job_result(job_id: str):
    job = query_one("SELECT * FROM analysis_jobs WHERE id = ?", (job_id,))
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")
    return {
        "job_id": job["id"],
        "status": job["status"],
        "result": json.loads(job["result_json"]) if job.get("result_json") else None,
    }