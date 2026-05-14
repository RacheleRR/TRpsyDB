"""
Expression router
GET /api/expression/{gene_symbol}   → GTEx v8 median expression per tissue, mapped to UBERON IDs

Note: GTEx v10 returns empty results for medianGeneExpression endpoint.
      GTEx v8 is stable and returns correct tissue expression data.
"""
import requests
from fastapi import APIRouter, HTTPException
from functools import lru_cache

router = APIRouter()

GTEX_API = "https://gtexportal.org/api/v2"

TISSUE_MAP = {
    "Adipose_Subcutaneous":                  {"uberon": "UBERON_0002190", "name": "Adipose (subcutaneous)",     "brain": False},
    "Adipose_Visceral_Omentum":              {"uberon": "UBERON_0010414", "name": "Adipose (visceral)",         "brain": False},
    "Adrenal_Gland":                         {"uberon": "UBERON_0002369", "name": "Adrenal gland",             "brain": False},
    "Artery_Aorta":                          {"uberon": "UBERON_0000947", "name": "Artery (aorta)",            "brain": False},
    "Artery_Coronary":                       {"uberon": "UBERON_0001621", "name": "Artery (coronary)",         "brain": False},
    "Artery_Tibial":                         {"uberon": "UBERON_0007610", "name": "Artery (tibial)",           "brain": False},
    "Bladder":                               {"uberon": "UBERON_0001255", "name": "Bladder",                   "brain": False},
    "Brain_Amygdala":                        {"uberon": "UBERON_0001876", "name": "Amygdala",                  "brain": True},
    "Brain_Anterior_cingulate_cortex_Ba24":  {"uberon": "UBERON_0001871", "name": "Anterior cingulate (BA24)","brain": True},
    "Brain_Caudate_basal_ganglia":           {"uberon": "UBERON_0001873", "name": "Caudate",                   "brain": True},
    "Brain_Cerebellar_Hemisphere":           {"uberon": "UBERON_0002245", "name": "Cerebellar hemisphere",     "brain": True},
    "Brain_Cerebellum":                      {"uberon": "UBERON_0002037", "name": "Cerebellum",                "brain": True},
    "Brain_Cortex":                          {"uberon": "UBERON_0000956", "name": "Cortex",                    "brain": True},
    "Brain_Frontal_Cortex_Ba9":              {"uberon": "UBERON_0001870", "name": "Frontal cortex (BA9)",      "brain": True},
    "Brain_Hippocampus":                     {"uberon": "UBERON_0001954", "name": "Hippocampus",               "brain": True},
    "Brain_Hypothalamus":                    {"uberon": "UBERON_0001898", "name": "Hypothalamus",              "brain": True},
    "Brain_Nucleus_accumbens_basal_ganglia": {"uberon": "UBERON_0001882", "name": "Nucleus accumbens",        "brain": True},
    "Brain_Putamen_basal_ganglia":           {"uberon": "UBERON_0001874", "name": "Putamen",                   "brain": True},
    "Brain_Spinal_cord_cervical_c_1":        {"uberon": "UBERON_0002240", "name": "Spinal cord (C1)",          "brain": True},
    "Brain_Substantia_nigra":               {"uberon": "UBERON_0002038", "name": "Substantia nigra",          "brain": True},
    "Breast_Mammary_Tissue":                 {"uberon": "UBERON_0001911", "name": "Breast",                    "brain": False},
    "Cells_Cultured_fibroblasts":            {"uberon": "UBERON_0015764", "name": "Fibroblasts",               "brain": False},
    "Colon_Sigmoid":                         {"uberon": "UBERON_0001159", "name": "Colon (sigmoid)",           "brain": False},
    "Colon_Transverse":                      {"uberon": "UBERON_0001153", "name": "Colon (transverse)",        "brain": False},
    "Esophagus_Gastroesophageal_Junction":   {"uberon": "UBERON_0001047", "name": "Esophagus (GEJ)",           "brain": False},
    "Esophagus_Mucosa":                      {"uberon": "UBERON_0001043", "name": "Esophagus (mucosa)",        "brain": False},
    "Esophagus_Muscularis":                  {"uberon": "UBERON_0001045", "name": "Esophagus (muscularis)",    "brain": False},
    "Heart_Atrial_Appendage":               {"uberon": "UBERON_0002079", "name": "Heart (atrial appendage)",  "brain": False},
    "Heart_Left_Ventricle":                  {"uberon": "UBERON_0002084", "name": "Heart (left ventricle)",    "brain": False},
    "Kidney_Cortex":                         {"uberon": "UBERON_0002113", "name": "Kidney cortex",             "brain": False},
    "Liver":                                 {"uberon": "UBERON_0002107", "name": "Liver",                     "brain": False},
    "Lung":                                  {"uberon": "UBERON_0002048", "name": "Lung",                      "brain": False},
    "Minor_Salivary_Gland":                  {"uberon": "UBERON_0001830", "name": "Salivary gland",            "brain": False},
    "Muscle_Skeletal":                       {"uberon": "UBERON_0001134", "name": "Skeletal muscle",           "brain": False},
    "Nerve_Tibial":                          {"uberon": "UBERON_0001321", "name": "Nerve (tibial)",            "brain": False},
    "Ovary":                                 {"uberon": "UBERON_0000992", "name": "Ovary",                     "brain": False},
    "Pancreas":                              {"uberon": "UBERON_0001264", "name": "Pancreas",                  "brain": False},
    "Pituitary":                             {"uberon": "UBERON_0000007", "name": "Pituitary",                 "brain": False},
    "Prostate":                              {"uberon": "UBERON_0002367", "name": "Prostate",                  "brain": False},
    "Skin_Not_Sun_Exposed_Suprapubic":       {"uberon": "UBERON_0036149", "name": "Skin (non-sun exposed)",    "brain": False},
    "Skin_Sun_Exposed_Lower_leg":            {"uberon": "UBERON_0013756", "name": "Skin (sun exposed)",        "brain": False},
    "Small_Intestine_Terminal_Ileum":        {"uberon": "UBERON_0002116", "name": "Small intestine",           "brain": False},
    "Spleen":                                {"uberon": "UBERON_0002106", "name": "Spleen",                    "brain": False},
    "Stomach":                               {"uberon": "UBERON_0000945", "name": "Stomach",                   "brain": False},
    "Testis":                                {"uberon": "UBERON_0000473", "name": "Testis",                    "brain": False},
    "Thyroid":                               {"uberon": "UBERON_0002046", "name": "Thyroid",                   "brain": False},
    "Uterus":                                {"uberon": "UBERON_0000995", "name": "Uterus",                    "brain": False},
    "Vagina":                                {"uberon": "UBERON_0000996", "name": "Vagina",                    "brain": False},
    "Whole_Blood":                           {"uberon": "UBERON_0000178", "name": "Whole blood",               "brain": False},
}


@lru_cache(maxsize=200)
def get_gencode_id(gene_symbol: str):
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
    except Exception as e:
        print(f"[GTEx] gene lookup failed for {gene_symbol}: {e}")
    return None


@lru_cache(maxsize=200)
def fetch_gtex_expression(gencode_id: str) -> list:
    try:
        r = requests.get(
            f"{GTEX_API}/expression/medianGeneExpression",
            params={"gencodeId": gencode_id, "datasetId": "gtex_v8", "itemsPerPage": 300},
            timeout=15,
        )
        r.raise_for_status()
        return r.json().get("data", [])
    except Exception as e:
        print(f"[GTEx] expression fetch failed for {gencode_id}: {e}")
        return []


@router.get("/{gene_symbol}")
def get_expression(gene_symbol: str):
    gene_symbol = gene_symbol.upper().strip()

    gencode_id = get_gencode_id(gene_symbol)
    if not gencode_id:
        raise HTTPException(status_code=404,
            detail=f"Gene '{gene_symbol}' not found in GTEx.")

    raw = fetch_gtex_expression(gencode_id)
    if not raw:
        raise HTTPException(status_code=404,
            detail=f"No expression data for '{gene_symbol}' in GTEx v8.")

    body, brain = [], []
    for entry in raw:
        tissue_id  = entry.get("tissueSiteDetailId", "")
        median_tpm = entry.get("median", 0.0)
        info       = TISSUE_MAP.get(tissue_id)
        if not info:
            continue
        record = {
            "tissue_gtex": tissue_id,
            "tissue_name": info["name"],
            "uberon_id":   info["uberon"],
            "median_tpm":  round(float(median_tpm), 3),
            "unit":        "TPM",
        }
        if info["brain"]:
            brain.append(record)
        else:
            body.append(record)

    body.sort(key=lambda x: -x["median_tpm"])
    brain.sort(key=lambda x: -x["median_tpm"])

    return {
        "gene_symbol":   gene_symbol,
        "gencode_id":    gencode_id,
        "dataset":       "GTEx v8",
        "n_tissues":     len(body) + len(brain),
        "max_body_tpm":  body[0]["median_tpm"]  if body  else 0,
        "max_brain_tpm": brain[0]["median_tpm"] if brain else 0,
        "body":          body,
        "brain":         brain,
        "top_body":      body[:8],
        "top_brain":     brain[:8],
    }
