"""
Step 8 — ingest TRxQTL files → tr_brain_qtl
Run from scripts/ folder:
python 08_ingest_trxqtl.py --input ../data/raw/trxqtl/ --db data/trs.db

Confirmed TRxQTL format (two formats found):

ROSMAP/NYGCALS/AnswerALS format:
  Col 1 (phenotype/gene): chr10:18531849:18533336:ENSG00000225527.1@RP11-383B4.4@-
  Col 2 (TR variant):     chr10_18558534_18558548_AATGG

GTEx format (same structure):
  Col 1: chr1:212858087:212858088:ENSG00000198468.7
  Col 2: chr1_212851126_212851131_AG

TR coordinates are in Col 2: chr_start_end_motif (underscore-separated)
Gene is parsed from Col 1.
"""
import sqlite3, gzip, csv, argparse, re
import pandas as pd
from pathlib import Path

def parse_filename(fname):
    """
    ROSMAP.DLPFC.TR-eQTL.txt.gz     → (ROSMAP, DLPFC, eQTL)
    NYGCALS.Cortex_Frontal.TR-sQTL  → (NYGCALS, Cortex_Frontal, sQTL)
    Brain_Amygdala.TR-eQTL.txt.gz   → (GTEx_v8, Brain_Amygdala, eQTL)
    """
    stem = re.sub(r'\.(txt\.gz|txt|gz)$', '', fname)
    qtl_match = re.search(r'TR-(\w+)', stem, re.IGNORECASE)
    qtl_type  = qtl_match.group(1) if qtl_match else "eQTL"
    stem_clean = re.sub(r'\.?TR-\w+', '', stem)

    if stem_clean.startswith("ROSMAP."):
        cohort = "ROSMAP"
        tissue = stem_clean.split(".", 2)[1] if "." in stem_clean else stem_clean
    elif stem_clean.startswith("NYGCALS."):
        cohort = "NYGCALS"
        tissue = stem_clean.split(".", 2)[1] if "." in stem_clean else stem_clean
    elif stem_clean.startswith("AnswerALS."):
        cohort = "AnswerALS"
        tissue = stem_clean.split(".", 2)[1] if "." in stem_clean else stem_clean
    else:
        cohort = "GTEx_v8"
        tissue = stem_clean

    return cohort, tissue, qtl_type

def parse_tr_from_variant(variant_str):
    """
    Parse TR coordinates from variant column.
    Format: chr10_18558534_18558548_AATGG → (chr10, 18558534, 18558548)
    """
    parts = str(variant_str).strip().split("_")
    if len(parts) < 3:
        return None, None, None
    chrom = parts[0]
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    try:
        start = int(parts[1])
        end   = int(parts[2])
        return chrom, start, end
    except:
        return None, None, None

def parse_gene_from_phenotype(pheno_str):
    """
    Parse gene from phenotype column.
    chr10:18531849:18533336:ENSG00000225527.1@RP11-383B4.4@- → ENSG00000225527.1
    chr1:212858087:212858088:ENSG00000198468.7 → ENSG00000198468.7
    """
    s = str(pheno_str).strip()
    # after last colon, before @
    parts = s.split(":")
    if len(parts) >= 4:
        gene_part = parts[3].split("@")[0]
        return gene_part
    return s

def safe_float(v):
    try:
        return float(v) if v and str(v).strip() not in ("","NA","nan",".") else None
    except:
        return None

def build_coord_index(db_path):
    """(chr, start, end) → tr_id for hg38 and hg19."""
    con = sqlite3.connect(db_path)
    rows = con.execute("""
        SELECT tr_id, chr, start, end, hg19_chr, hg19_start, hg19_end
        FROM trs
    """).fetchall()
    con.close()

    index = {}
    for tr_id, c38, s38, e38, c19, s19, e19 in rows:
        index[(str(c38), int(s38), int(e38))]     = tr_id
        raw = str(c38).replace("chr","")
        index[(raw, int(s38), int(e38))]           = tr_id
        if c19 and s19 and e19:
            index[(str(c19), int(s19), int(e19))]  = tr_id
            raw19 = str(c19).replace("chr","")
            index[(raw19, int(s19), int(e19))]     = tr_id
    print(f"  Coord index: {len(index):,} entries")
    return index

def ingest_file(fpath, cohort, tissue, qtl_type, coord_index):
    rows = []
    opener = gzip.open if str(fpath).endswith(".gz") else open
    try:
        with opener(fpath, "rt", errors="replace") as f:
            # peek at header
            first = f.readline().strip()
            f.seek(0)

            # detect columns from header
            cols = first.split("\t")
            # TRxQTL files have no header — first line is data
            # detect: if first token looks like chr coordinate → no header
            has_header = not (cols[0].startswith("chr") or
                              re.match(r'chr\w+:\d+', cols[0]))

            if has_header:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    # try standard column names
                    pheno   = (row.get("phenotype_id") or row.get("gene_id") or
                               row.get("phenotype") or list(row.values())[0])
                    variant = (row.get("variant_id") or row.get("tr_id") or
                               list(row.values())[1])
                    p_val   = safe_float(row.get("pval_nominal") or
                                         row.get("p_value") or row.get("p_nom"))
                    beta    = safe_float(row.get("slope") or row.get("beta"))
                    se      = safe_float(row.get("slope_se") or row.get("se"))
                    q_val   = safe_float(row.get("qval") or row.get("q_value"))
                    gene_id = parse_gene_from_phenotype(pheno)
                    chrom, start, end = parse_tr_from_variant(variant)
                    if not chrom: continue
                    tr_id = coord_index.get((chrom, start, end))
                    if not tr_id: continue
                    rows.append((tr_id, cohort, tissue, qtl_type,
                                 gene_id, beta, se, p_val, q_val))
            else:
                # no header — positional columns
                # col0=phenotype, col1=variant, col2=p_value, col3=beta
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 2: continue
                    pheno   = parts[0]
                    variant = parts[1]
                    p_val   = safe_float(parts[2]) if len(parts) > 2 else None
                    beta    = safe_float(parts[3]) if len(parts) > 3 else None
                    gene_id = parse_gene_from_phenotype(pheno)
                    chrom, start, end = parse_tr_from_variant(variant)
                    if not chrom: continue
                    tr_id = coord_index.get((chrom, start, end))
                    if not tr_id: continue
                    rows.append((tr_id, cohort, tissue, qtl_type,
                                 gene_id, beta, None, p_val, None))
    except Exception as e:
        print(f"  ERROR: {e}")
    return rows

def ingest(input_dir, db_path):
    print("Building coordinate index ...")
    coord_index = build_coord_index(db_path)

    files = sorted(Path(input_dir).glob("*.gz")) + \
            sorted(Path(input_dir).glob("*.txt"))
    print(f"Found {len(files)} TRxQTL files\n")

    con = sqlite3.connect(db_path)
    con.execute("PRAGMA journal_mode = WAL")
    con.execute("PRAGMA synchronous = NORMAL")
    cur = con.cursor()

    total = 0
    for fpath in files:
        cohort, tissue, qtl_type = parse_filename(fpath.name)
        print(f"  {fpath.name}")
        print(f"    → {cohort} / {tissue} / {qtl_type}", end=" ... ")

        rows = ingest_file(fpath, cohort, tissue, qtl_type, coord_index)
        if rows:
            cur.executemany("""
                INSERT INTO tr_brain_qtl
                (tr_id, cohort, tissue, qtl_type, gene_id,
                 beta, se, p_value, q_value)
                VALUES (?,?,?,?,?,?,?,?,?)
            """, rows)
            con.commit()
            total += len(rows)
            print(f"{len(rows):,} rows matched")
        else:
            print("0 matched")

    n = cur.execute("SELECT COUNT(*) FROM tr_brain_qtl").fetchone()[0]
    print(f"\nDone. tr_brain_qtl: {n:,} rows total")

    # summary by cohort
    summary = cur.execute("""
        SELECT cohort, qtl_type, COUNT(*) as n
        FROM tr_brain_qtl
        GROUP BY cohort, qtl_type
        ORDER BY cohort, qtl_type
    """).fetchall()
    for cohort, qt, n in summary:
        print(f"  {cohort:15s} {qt:8s} {n:,}")
    con.close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--db",    default="data/trs.db")
    args = ap.parse_args()
    ingest(args.input, args.db)