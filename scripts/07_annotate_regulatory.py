"""
Step 7 — annotate regulatory features for all TRs
Run from scripts/ folder:
python 07_annotate_regulatory.py \
    --db        data/trs.db \
    --epd       /home/rachele/TRpsyDB/data/raw/regulatory/Hs_EPDnew.bed \
    --ccre      /home/rachele/TRpsyDB/data/raw/regulatory/cCREs_hg38.bed \
    --sedb_se   /home/rachele/TRpsyDB/data/raw/regulatory/SE.bed \
    --sedb_te   /home/rachele/TRpsyDB/data/raw/regulatory/SE_te.bed \
    --sedb_meta /home/rachele/TRpsyDB/data/raw/regulatory/Human_sample_information_sedb3.txt \
    --refflat   /home/rachele/TRpsyDB/data/raw/refflat/refFlat.txt \
    --introns   /home/rachele/TRpsyDB/data/raw/refflat/hg38_intron_refFlat.txt \
    --exons     /home/rachele/TRpsyDB/data/raw/refflat/hg38_exon_refFlat.txt \
    --appris    /home/rachele/TRpsyDB/data/raw/refflat/appris_data.principal.hg38.txt
"""
import sqlite3, argparse, subprocess
import pandas as pd
import pyranges as pr
from pathlib import Path

# ── helpers ───────────────────────────────────────────────────────────────────
def load_trs_df(db_path):
    con = sqlite3.connect(db_path)
    df = pd.read_sql("SELECT tr_id, chr, start, end FROM trs", con)
    con.close()
    return df

def trs_to_pr(df):
    tmp = df.rename(columns={"chr":"Chromosome","start":"Start","end":"End"})
    tmp["Chromosome"] = tmp["Chromosome"].astype(str)
    return pr.PyRanges(tmp)

def pr_join_back(trs_df, right_pr, right_cols):
    """
    Inner join trs_pr with right_pr, then left-merge result back to
    full trs_df so every TR has a row (NaN for non-overlapping).
    Returns DataFrame with tr_id + right_cols.
    """
    trs_pr = trs_to_pr(trs_df)
    # inner join — only overlapping TRs
    joined = trs_pr.join(right_pr)
    if len(joined) == 0:
        result = trs_df[["tr_id"]].copy()
        for c in right_cols:
            result[c] = None
        return result
    hit_df = joined.df[["tr_id"] + right_cols].copy()
    # keep best hit per TR (first non-null)
    hit_df = hit_df.groupby("tr_id").first().reset_index()
    # merge back to full list
    result = trs_df[["tr_id"]].merge(hit_df, on="tr_id", how="left")
    return result

# ── EPD promoters ─────────────────────────────────────────────────────────────
def annotate_epd(trs_df, epd_path):
    print("Annotating EPD promoters ...")
    epd = pd.read_csv(epd_path, sep=r'\s+', header=None, low_memory=False,
        names=["Chromosome","Start","End","Name","Score","Strand","ThickStart","ThickEnd"])
    epd["Chromosome"] = epd["Chromosome"].astype(str)
    epd_pr = pr.PyRanges(epd[["Chromosome","Start","End","Name"]])

    result = pr_join_back(trs_df, epd_pr, ["Name"])
    result = result.rename(columns={"Name":"promoter_id"})
    result["in_promoter"] = result["promoter_id"].notna().astype(int)
    print(f"  {result['in_promoter'].sum():,} / {len(result):,} TRs overlap a promoter")
    return result

# ── ENCODE cCREs ──────────────────────────────────────────────────────────────
# Format: chr start end accession_d accession_e  class[,class]
def annotate_ccre(trs_df, ccre_path):
    print("Annotating ENCODE cCREs ...")
    ccre = pd.read_csv(ccre_path, sep="\t", header=None, low_memory=False,
        names=["Chromosome","Start","End","ccre_d","ccre_id","ccre_raw"])
    ccre["Chromosome"] = ccre["Chromosome"].astype(str)
    ccre["ccre_class"] = ccre["ccre_raw"].str.split(",").str[0]
    ccre_pr = pr.PyRanges(ccre[["Chromosome","Start","End","ccre_id","ccre_class"]])

    result = pr_join_back(trs_df, ccre_pr, ["ccre_id","ccre_class"])
    result["in_any_ccre"]   = result["ccre_class"].notna().astype(int)
    result["in_brain_ccre"] = 0  # TODO v2: add brain biosample filter
    print(f"  {result['in_any_ccre'].sum():,} / {len(result):,} TRs overlap a cCRE")
    return result

# ── SEdb brain cell_ids ───────────────────────────────────────────────────────
BRAIN_KEYWORDS = (
    "brain|neuro|cortex|cerebel|hippoc|amygdala|striatum|"
    "substantia|caudate|putamen|nucleus accumbens|frontal|temporal|"
    "occipital|parietal|thalamus|hypothalamus|midbrain|pons|medulla"
)

def load_brain_cells(meta_path):
    if not meta_path or not Path(meta_path).exists():
        return set()
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    # confirmed cols: Sample ID | Tissue type
    if "Tissue type" in meta.columns and "Sample ID" in meta.columns:
        brain = meta[meta["Tissue type"].str.contains(
            BRAIN_KEYWORDS, case=False, na=False)]
        ids = set(brain["Sample ID"].str.strip())
        print(f"  Brain cell_ids from metadata: {len(ids)}")
        return ids
    return set()

# ── SEdb super enhancers ──────────────────────────────────────────────────────
# Confirmed format: quoted TSV
# "cell_id" "se_id" "se_chr" "se_start" "se_end" ... "se_gene_overlap" "se_gene_closest" "se_gene_ABC"
def annotate_sedb_se(trs_df, se_path, brain_cells):
    print("Annotating SEdb super enhancers ...")
    se = pd.read_csv(se_path, sep="\t", low_memory=False, quotechar='"', dtype=str)
    se.columns = [c.strip('"') for c in se.columns]
    for c in se.columns:
        se[c] = se[c].str.strip('"').str.strip()

    se = se.rename(columns={"se_chr":"Chromosome","se_start":"Start","se_end":"End"})
    se["Start"] = pd.to_numeric(se["Start"], errors="coerce")
    se["End"]   = pd.to_numeric(se["End"],   errors="coerce")
    se = se.dropna(subset=["Start","End"])
    se["Chromosome"] = se["Chromosome"].astype(str)
    se["is_brain"]   = se["cell_id"].isin(brain_cells).astype(int)

    # gene columns — take first value per SE
    gene_cols = [c for c in ["se_gene_overlap","se_gene_closest","se_gene_ABC","se_cas_value"]
                 if c in se.columns]
    keep = ["Chromosome","Start","End","cell_id","is_brain"] + gene_cols
    se_pr = pr.PyRanges(se[keep])

    joined = trs_to_pr(trs_df).join(se_pr)
    if len(joined) == 0:
        result = trs_df[["tr_id"]].copy()
        for c in ["in_any_super_enhancer","in_brain_super_enhancer",
                  "n_super_enhancer_tissues","super_enhancer_tissues",
                  "SE_gene_overlap","SE_gene_closest","SE_gene_ABC","SE_cas_value"]:
            result[c] = None
        return result

    df = joined.df.copy()
    agg = df.groupby("tr_id").agg(
        in_any_super_enhancer   =("cell_id",  lambda x: int(x.notna().any())),
        in_brain_super_enhancer =("is_brain", lambda x: int((x.astype(int)==1).any())),
        n_super_enhancer_tissues=("cell_id",  "nunique"),
        super_enhancer_tissues  =("cell_id",  lambda x: ",".join(x.dropna().unique()[:5])),
    )
    for g in gene_cols:
        agg[g.replace("se_","SE_")] = df.groupby("tr_id")[g].first()
    agg = agg.reset_index()

    result = trs_df[["tr_id"]].merge(agg, on="tr_id", how="left")
    for c in ["in_any_super_enhancer","in_brain_super_enhancer","n_super_enhancer_tissues"]:
        result[c] = result[c].fillna(0).astype(int)
    print(f"  {result['in_any_super_enhancer'].sum():,} TRs in any SE")
    print(f"  {result['in_brain_super_enhancer'].sum():,} TRs in brain SE")
    return result

# ── SEdb typical enhancers ────────────────────────────────────────────────────
# Confirmed format: no quotes
# te_chr te_start te_end cell_id te_id te_rank te_ele_num te_CONSTITUENT_SIZE te_cas_value te_con_value isSUPER
def annotate_sedb_te(trs_df, te_path, brain_cells):
    print("Annotating SEdb typical enhancers ...")
    te = pd.read_csv(te_path, sep="\t", low_memory=False, dtype=str)
    te = te.rename(columns={"te_chr":"Chromosome","te_start":"Start","te_end":"End"})
    te["Start"] = pd.to_numeric(te["Start"], errors="coerce")
    te["End"]   = pd.to_numeric(te["End"],   errors="coerce")
    te = te.dropna(subset=["Start","End"])
    te["Chromosome"] = te["Chromosome"].astype(str)
    te["is_brain"]   = te["cell_id"].isin(brain_cells).astype(int)

    te_pr = pr.PyRanges(te[["Chromosome","Start","End","cell_id","is_brain"]])
    joined = trs_to_pr(trs_df).join(te_pr)
    if len(joined) == 0:
        result = trs_df[["tr_id"]].copy()
        for c in ["in_any_enhancer","in_brain_enhancer",
                  "n_enhancer_tissues","enhancer_tissues"]:
            result[c] = None
        return result

    df = joined.df.copy()
    agg = df.groupby("tr_id").agg(
        in_any_enhancer   =("cell_id",  lambda x: int(x.notna().any())),
        in_brain_enhancer =("is_brain", lambda x: int((x.astype(int)==1).any())),
        n_enhancer_tissues=("cell_id",  "nunique"),
        enhancer_tissues  =("cell_id",  lambda x: ",".join(x.dropna().unique()[:5])),
    ).reset_index()

    result = trs_df[["tr_id"]].merge(agg, on="tr_id", how="left")
    for c in ["in_any_enhancer","in_brain_enhancer","n_enhancer_tissues"]:
        result[c] = result[c].fillna(0).astype(int)
    print(f"  {result['in_any_enhancer'].sum():,} TRs in any TE")
    print(f"  {result['in_brain_enhancer'].sum():,} TRs in brain TE")
    return result

# ── TSS + SJ distances via R ──────────────────────────────────────────────────
def compute_tss_sj(db_path, refflat, introns, appris108=None, appris109=None):
    print("Computing TSS and SJ distances via R ...")
    Path("data/processed").mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(db_path)
    df  = pd.read_sql("SELECT tr_id, chr, start, end FROM trs", con)
    con.close()
    trs_bed = "data/processed/trs_for_tss.bed"
    df.to_csv(trs_bed, sep="\t", index=False, header=False)

    # NOTE: introns are RefSeq-based (NM_... isoform IDs)
    #       APPRIS uses Ensembl IDs (ENST...) — incompatible filter, skip it
    #       refFlat already represents canonical transcripts well enough
    r_script = r"""
library(dplyr)
library(GenomicRanges)
options(warn=1)

args        <- commandArgs(trailingOnly=TRUE)
trs_bed     <- args[1]
refflat_f   <- args[2]
introns_f   <- args[3]
appris108_f <- args[4]
appris109_f <- args[5]
out_tss     <- args[6]
out_sj      <- args[7]

cat("Loading trs ...\n")
trs <- read.delim(trs_bed, header=FALSE,
         col.names=c("tr_id","chrom","start","end"))
cat("  ", nrow(trs), "TRs\n")

cat("Loading refFlat ...\n")
refflat <- read.delim(refflat_f, header=FALSE,
    col.names=c("geneName","name","chrom","strand","txStart","txEnd",
                "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))
cat("  ", nrow(refflat), "transcripts\n")

cat("Loading introns ...\n")
introns <- read.delim(introns_f, header=TRUE, stringsAsFactors=FALSE)
cat("  ", nrow(introns), "introns before APPRIS filter\n")

# load APPRIS RefSeq files (no header, col3 = transcript ID)
appris_list <- list()
if (nchar(appris108_f) > 0 && file.exists(appris108_f)) {
    appris_list[[1]] <- read.delim(appris108_f, header=FALSE, stringsAsFactors=FALSE)
    cat("  Loaded appris108:", nrow(appris_list[[1]]), "\n")
}
if (nchar(appris109_f) > 0 && file.exists(appris109_f)) {
    appris_list[[2]] <- read.delim(appris109_f, header=FALSE, stringsAsFactors=FALSE)
    cat("  Loaded appris109:", nrow(appris_list[[2]]), "\n")
}
if (length(appris_list) > 0) {
    appris <- do.call(rbind, appris_list)
    # V3 = transcript ID e.g. NM_005465.4 — strip version number
    appris_ids <- sapply(strsplit(appris$V3, "\\."), "[", 1)
    # filter to PRINCIPAL isoforms only
    is_principal <- grepl("PRINCIPAL", appris$V5)
    appris_ids <- appris_ids[is_principal]
    introns <- introns[introns$isoform %in% appris_ids, ]
    cat("  Introns after APPRIS filter:", nrow(introns), "\n")
} else {
    cat("  No APPRIS files found — using all introns\n")
}

# TSS distance
tss_df <- refflat %>%
    mutate(tss = ifelse(strand=="+", txStart, txEnd)) %>%
    select(geneName, chrom, tss) %>%
    distinct()

trs_gr <- GRanges(trs$chrom, IRanges(trs$start + 1, trs$end))
tss_gr <- GRanges(tss_df$chrom, IRanges(tss_df$tss, tss_df$tss),
                  gene=tss_df$geneName)

cat("Computing TSS distances ...\n")
hits <- distanceToNearest(trs_gr, tss_gr)
tss_out <- data.frame(
    tr_id            = trs$tr_id[queryHits(hits)],
    dist_tss_bp      = mcols(hits)$distance,
    nearest_tss_gene = tss_gr$gene[subjectHits(hits)]
)
write.table(tss_out, out_tss, sep="\t", row.names=FALSE, quote=FALSE)
cat("TSS done:", nrow(tss_out), "\n")

# SJ distance — distance to nearest splice junction (intron boundaries)
cat("Computing SJ distances ...\n")
sj_pos <- c(as.integer(introns$start), as.integer(introns$end))
sj_chr <- c(as.character(introns$chrom), as.character(introns$chrom))
valid  <- !is.na(sj_pos) & nchar(sj_chr) > 0
sj_gr  <- GRanges(sj_chr[valid], IRanges(sj_pos[valid], sj_pos[valid]))

cat("  SJ positions:", length(sj_gr), "\n")
hits_sj <- distanceToNearest(trs_gr, sj_gr)
sj_out  <- data.frame(
    tr_id          = trs$tr_id[queryHits(hits_sj)],
    dist_splice_bp = mcols(hits_sj)$distance
)
write.table(sj_out, out_sj, sep="\t", row.names=FALSE, quote=FALSE)
cat("SJ done:", nrow(sj_out), "\n")
"""
    rfile = "data/processed/07b_tss_sj.R"
    with open(rfile, "w") as f:
        f.write(r_script)

    out_tss = "data/processed/tss_distances.tsv"
    out_sj  = "data/processed/sj_distances.tsv"

    # remove old output files so we can detect failure
    for p in [out_tss, out_sj]:
        Path(p).unlink(missing_ok=True)

    res = subprocess.run(
        ["Rscript", rfile, trs_bed, refflat, introns,
         appris108 or "", appris109 or "", out_tss, out_sj],
        capture_output=True, text=True
    )
    print(res.stdout)
    if res.returncode != 0:
        print("R STDERR:", res.stderr[-3000:])
        return None, None

    tss_df = pd.read_csv(out_tss, sep="\t") if Path(out_tss).exists() else None
    sj_df  = pd.read_csv(out_sj,  sep="\t") if Path(out_sj).exists() else None
    if tss_df is not None: print(f"  TSS rows: {len(tss_df):,}")
    if sj_df  is not None: print(f"  SJ rows:  {len(sj_df):,}")
    return tss_df, sj_df

# ── Write to tr_regulatory ────────────────────────────────────────────────────
def write_regulatory(db_path, epd_df, ccre_df, se_df, te_df, tss_df, sj_df):
    print("Writing to tr_regulatory ...")
    con = sqlite3.connect(db_path)

    # recreate tr_regulatory with correct schema
    con.execute("DROP TABLE IF EXISTS tr_regulatory")
    con.execute("""
        CREATE TABLE tr_regulatory (
            tr_id                    TEXT PRIMARY KEY REFERENCES trs(tr_id),
            dist_tss_bp              INTEGER,
            nearest_tss_gene         TEXT,
            dist_splice_bp           INTEGER,
            in_promoter              INTEGER DEFAULT 0,
            promoter_id              TEXT,
            ccre_class               TEXT,
            ccre_id                  TEXT,
            in_any_ccre              INTEGER DEFAULT 0,
            in_brain_ccre            INTEGER DEFAULT 0,
            in_any_enhancer          INTEGER DEFAULT 0,
            in_brain_enhancer        INTEGER DEFAULT 0,
            enhancer_tissues         TEXT,
            n_enhancer_tissues       INTEGER DEFAULT 0,
            in_any_super_enhancer    INTEGER DEFAULT 0,
            in_brain_super_enhancer  INTEGER DEFAULT 0,
            super_enhancer_tissues   TEXT,
            n_super_enhancer_tissues INTEGER DEFAULT 0,
            SE_gene_overlap          TEXT,
            SE_gene_closest          TEXT,
            SE_gene_ABC              TEXT,
            SE_cas_value             REAL
        )
    """)
    con.commit()

    # start with full trs list
    base = pd.read_sql("SELECT tr_id FROM trs", con)
    result = base.copy()

    for df in [epd_df, ccre_df, se_df, te_df, tss_df, sj_df]:
        if df is not None:
            result = result.merge(df, on="tr_id", how="left")

    # fill boolean NAs
    bool_cols = ["in_promoter","in_any_ccre","in_brain_ccre",
                 "in_any_enhancer","in_brain_enhancer",
                 "in_any_super_enhancer","in_brain_super_enhancer"]
    for c in bool_cols:
        if c in result.columns:
            result[c] = result[c].fillna(0).astype(int)

    # only keep columns that exist in the table
    cur = con.cursor()
    table_cols = [r[1] for r in cur.execute("PRAGMA table_info(tr_regulatory)")]
    result = result[[c for c in result.columns if c in table_cols]]

    chunk = 50000
    for s in range(0, len(result), chunk):
        result.iloc[s:s+chunk].to_sql(
            "tr_regulatory", con, if_exists="append", index=False)
        con.commit()
        print(f"  {min(s+chunk, len(result)):,} / {len(result):,}", end="\r")
    print()

    n = con.execute("SELECT COUNT(*) FROM tr_regulatory").fetchone()[0]
    n_prom = con.execute("SELECT COUNT(*) FROM tr_regulatory WHERE in_promoter=1").fetchone()[0]
    n_ccre = con.execute("SELECT COUNT(*) FROM tr_regulatory WHERE in_any_ccre=1").fetchone()[0]
    n_se   = con.execute("SELECT COUNT(*) FROM tr_regulatory WHERE in_any_super_enhancer=1").fetchone()[0]
    n_te   = con.execute("SELECT COUNT(*) FROM tr_regulatory WHERE in_any_enhancer=1").fetchone()[0]
    print(f"\nDone. tr_regulatory: {n:,} rows")
    print(f"  in_promoter:             {n_prom:,}")
    print(f"  in_any_ccre:             {n_ccre:,}")
    print(f"  in_any_super_enhancer:   {n_se:,}")
    print(f"  in_any_enhancer:         {n_te:,}")
    con.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db");       ap.add_argument("--epd")
    ap.add_argument("--ccre");     ap.add_argument("--sedb_se", default=None)
    ap.add_argument("--sedb_te",   default=None)
    ap.add_argument("--sedb_meta", default=None)
    ap.add_argument("--refflat");  ap.add_argument("--introns")
    ap.add_argument("--exons")
    ap.add_argument("--appris108", default=None,
                    help="appris_data.principal.refseq108.hg38.txt")
    ap.add_argument("--appris109", default=None,
                    help="appris_data.principal.refseq109.hg38.txt")
    args = ap.parse_args()

    trs_df = load_trs_df(args.db)
    print(f"Loaded {len(trs_df):,} TRs from DB")

    brain_cells = load_brain_cells(args.sedb_meta)

    epd_df  = annotate_epd(trs_df, args.epd)
    ccre_df = annotate_ccre(trs_df, args.ccre)
    se_df   = annotate_sedb_se(trs_df, args.sedb_se, brain_cells) if args.sedb_se else None
    te_df   = annotate_sedb_te(trs_df, args.sedb_te, brain_cells) if args.sedb_te else None

    tss_df, sj_df = compute_tss_sj(
        args.db, args.refflat, args.introns,
        args.appris108, args.appris109
    )
    write_regulatory(args.db, epd_df, ccre_df, se_df, te_df, tss_df, sj_df)

if __name__ == "__main__":
    main()