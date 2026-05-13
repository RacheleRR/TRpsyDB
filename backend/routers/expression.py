"""
Expression router
GET /api/expression/{gene_symbol}   → GTEx v10 expression per tissue, mapped to UBERON IDs
"""
import requests
from fastapi import APIRouter, HTTPException
from functools import lru_cache

router = APIRouter()

GTEX_API = "https://gtexportal.org/api/v2"

# ── GTEx tissue name → UBERON ID mapping ─────────────────────
GTEX_TO_UBERON = {
    "Brain"                                      : "UBERON_0000955",
    "Lung"                                        : "UBERON_0002048",
    "Heart - Left Ventricle"                      : "UBERON_0002084",
    "Heart - Atrial Appendage"                    : "UBERON_0002079",
    "Liver"                                       : "UBERON_0002107",
    "Kidney - Cortex"                             : "UBERON_0002113",
    "Kidney - Medulla"                            : "UBERON_0000014",
    "Colon - Sigmoid"                             : "UBERON_0001155",
    "Colon - Transverse"                          : "UBERON_0001153",
    "Muscle - Skeletal"                           : "UBERON_0001134",
    "Skin - Sun Exposed (Lower leg)"              : "UBERON_0001013",
    "Adipose - Subcutaneous"                      : "UBERON_0002190",
    "Adipose - Visceral (Omentum)"                : "UBERON_0010414",
    "Thyroid"                                     : "UBERON_0002046",
    "Testis"                                      : "UBERON_0000473",
    "Whole Blood"                                 : "UBERON_0000178",
    "Stomach"                                     : "UBERON_0000945",
    "Spleen"                                      : "UBERON_0002106",
    "Pancreas"                                    : "UBERON_0001264",
    "Adrenal Gland"                               : "UBERON_0002369",
    "Pituitary"                                   : "UBERON_0000007",
    "Small Intestine - Terminal Ileum"            : "UBERON_0002116",
    "Esophagus - Mucosa"                          : "UBERON_0001043",
    "Esophagus - Muscularis"                      : "UBERON_0001045",
    "Esophagus - Gastroesophageal Junction"       : "UBERON_0001047",
    "Bladder"                                     : "UBERON_0001255",
    "Nerve - Tibial"                              : "UBERON_0001021",
    "Ovary"                                       : "UBERON_0000992",
    "Uterus"                                      : "UBERON_0000995",
    "Prostate"                                    : "UBERON_0002367",
    "Breast - Mammary Tissue"                     : "UBERON_0001911",
    "Vagina"                                      : "UBERON_0000996",
    "Artery - Aorta"                              : "UBERON_0000947",
    "Artery - Coronary"                           : "UBERON_0001621",
    "Artery - Tibial"                             : "UBERON_0007610",
    "Minor Salivary Gland"                        : "UBERON_0001830",
    # Brain regions
    "Brain - Frontal Cortex (BA9)"                : "UBERON_0001870",
    "Brain - Anterior cingulate cortex (BA24)"    : "UBERON_0003027",
    "Brain - Caudate (basal ganglia)"             : "UBERON_0001873",
    "Brain - Cerebellar Hemisphere"               : "UBERON_0002245",
    "Brain - Cerebellum"                          : "UBERON_0002285",
    "Brain - Cortex"                              : "UBERON_0000956",
    "Brain - Hippocampus"                         : "UBERON_0002421",
    "Brain - Hypothalamus"                        : "UBERON_0001898",
    "Brain - Nucleus accumbens (basal ganglia)"   : "UBERON_0001875",
    "Brain - Putamen (basal ganglia)"             : "UBERON_0001874",
    "Brain - Spinal cord (cervical c-1)"          : "UBERON_0002360",
    "Brain - Substantia nigra"                    : "UBERON_0002038",
    "Brain - Amygdala"                            : "UBERON_0001876",
}

BRAIN_TISSUES = {k for k in GTEX_TO_UBERON if k.startswith("Brain")}


@lru_cache(maxsize=200)
def get_ensembl_id(gene_symbol: str) -> str | None:
    """Resolve gene symbol → Ensembl versioned ID via GTEx API."""
    try:
        r = requests.get(
            f"{GTEX_API}/reference/gene",
            params={"geneId": gene_symbol, "gencodeVersion": "v26",
                    "genomeBuild": "GRCh38/hg38", "pageSize": 5},
            timeout=10,
        )
        r.raise_for_status()
        data = r.json().get("data", [])
        if data:
            return data[0].get("gencodeId")
    except Exception:
        pass
    return None


@lru_cache(maxsize=200)
def fetch_gtex_expression(gencode_id: str) -> list:
    """Fetch median TPM per tissue from GTEx v10."""
    try:
        r = requests.get(
            f"{GTEX_API}/expression/geneExpression",
            params={"gencodeId": gencode_id, "datasetId": "gtex_v10"},
            timeout=15,
        )
        r.raise_for_status()
        return r.json().get("data", [])
    except Exception:
        return []


@router.get("/{gene_symbol}")
def get_expression(gene_symbol: str):
    """
    Returns GTEx v10 median TPM per tissue, mapped to UBERON IDs.
    Separates body tissues and brain regions for frontend rendering.
    """
    gene_symbol = gene_symbol.upper().strip()

    gencode_id = get_ensembl_id(gene_symbol)
    if not gencode_id:
        raise HTTPException(
            status_code=404,
            detail=f"Gene '{gene_symbol}' not found in GTEx. Check the gene symbol."
        )

    raw = fetch_gtex_expression(gencode_id)
    if not raw:
        raise HTTPException(
            status_code=404,
            detail=f"No expression data found for '{gene_symbol}' in GTEx v10."
        )

    body_data   = {}   # uberon_id → {tissue_name, median_tpm, unit}
    brain_data  = {}

    for entry in raw:
        tissue_name = entry.get("tissueSiteDetailId", "").replace("_", " - ").replace("  ", " ")
        # GTEx uses underscore-separated tissue IDs; try to match to our mapping
        median_tpm  = entry.get("median", 0.0)

        # Try direct match first, then normalised match
        uberon = GTEX_TO_UBERON.get(tissue_name)
        if not uberon:
            # GTEx API returns IDs like "Brain_Frontal_Cortex_Ba9"
            # try to find a close match
            for gtex_name, uid in GTEX_TO_UBERON.items():
                if gtex_name.lower().replace(" ", "_") == tissue_name.lower().replace(" ", "_"):
                    uberon = uid
                    tissue_name = gtex_name
                    break

        if not uberon:
            continue

        record = {
            "tissue_name": tissue_name,
            "uberon_id":   uberon,
            "median_tpm":  round(median_tpm, 3),
            "unit":        "TPM",
        }

        if tissue_name in BRAIN_TISSUES:
            brain_data[uberon] = record
        else:
            body_data[uberon] = record

    # Compute max TPM for normalisation (frontend uses this for colour scale)
    all_tpms     = [v["median_tpm"] for v in {**body_data, **brain_data}.values()]
    max_body_tpm = max((v["median_tpm"] for v in body_data.values()), default=0)
    max_brain_tpm= max((v["median_tpm"] for v in brain_data.values()), default=0)

    # Top expressed tissues (for the ranked list)
    top_body  = sorted(body_data.values(),  key=lambda x: -x["median_tpm"])[:8]
    top_brain = sorted(brain_data.values(), key=lambda x: -x["median_tpm"])[:8]

    return {
        "gene_symbol":   gene_symbol,
        "gencode_id":    gencode_id,
        "dataset":       "GTEx v10",
        "n_tissues":     len(body_data) + len(brain_data),
        "max_body_tpm":  round(max_body_tpm, 3),
        "max_brain_tpm": round(max_brain_tpm, 3),
        "body":          body_data,
        "brain":         brain_data,
        "top_body":      top_body,
        "top_brain":     top_brain,
    }
