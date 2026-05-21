"""
Step 6 — ingest cleaned RExPRT output → tr_rexprt table
Run: python 06_ingest_rexprt.py --input data/processed/rexprt_cleaned.tsv --db data/trs.db
"""
import sqlite3, csv, argparse
from pathlib import Path

def safe_float(v):
    try: return float(v) if v and str(v).strip() not in ("", "NA", "NaN", ".") else None
    except: return None

def safe_int(v):
    try: return int(float(v)) if v and str(v).strip() not in ("", "NA", "NaN", ".") else None
    except: return None

def ingest(input_path, db_path):
    print(f"Reading {input_path} ...")
    with open(input_path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"  {len(rows):,} rows")

    con = sqlite3.connect(db_path)
    con.execute("PRAGMA foreign_keys = ON")
    con.execute("PRAGMA journal_mode = WAL")
    cur = con.cursor()

    # get all valid tr_ids for checking
    valid_ids = {r[0] for r in cur.execute("SELECT tr_id FROM trs").fetchall()}

    inserted = skipped = 0
    rexprt_rows = []

    for row in rows:
        tr_id = row.get("tr_id", "").strip()
        if not tr_id or tr_id not in valid_ids:
            skipped += 1
            continue

        rexprt_rows.append((
            tr_id,
            safe_float(row.get("ensembleScore")),
            safe_int(row.get("ensembleBinary")),
            safe_float(row.get("ensembleMax")),
            safe_float(row.get("SVM")),
            safe_float(row.get("XGB")),
            safe_float(row.get("pLi")),
            safe_float(row.get("loeuf")),
            safe_float(row.get("gc_content")),
            safe_float(row.get("per_g")),
            safe_float(row.get("per_c")),
            safe_float(row.get("per_a")),
            safe_float(row.get("per_t")),
            safe_int(row.get("eSTR")),
            safe_int(row.get("opReg")),
            safe_int(row.get("promoter")),
            safe_int(row.get("UTR_5")),
            safe_int(row.get("UTR_3")),
            row.get("HG19_ID", ""),
        ))
        inserted += 1

    cur.executemany("""
        INSERT OR REPLACE INTO tr_rexprt (
            tr_id, ensembleScore, ensembleBinary, ensembleMax,
            SVM, XGB, pLi, loeuf,
            gc_content, per_g, per_c, per_a, per_t,
            eSTR, opReg, promoter, UTR_5, UTR_3,
            HG19_ID
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, rexprt_rows)

    # also update hg19 coords in trs from HG19_ID
    for row in rows:
        tr_id  = row.get("tr_id", "").strip()
        hg19   = row.get("HG19_ID", "")
        if tr_id and hg19 and "#" in hg19:
            parts = hg19.split("#")
            if len(parts) == 3:
                cur.execute("""
                    UPDATE trs SET hg19_chr=?, hg19_start=?, hg19_end=?
                    WHERE tr_id=? AND hg19_chr IS NULL
                """, (parts[0], safe_int(parts[1]), safe_int(parts[2]), tr_id))

    con.commit()

    n_rexprt = cur.execute("SELECT COUNT(*) FROM tr_rexprt").fetchone()[0]
    n_path   = cur.execute(
        "SELECT COUNT(*) FROM tr_rexprt WHERE ensembleBinary = 1"
    ).fetchone()[0]
    print(f"\nDone.")
    print(f"  Inserted into tr_rexprt: {inserted:,}")
    print(f"  Skipped (no tr_id):      {skipped}")
    print(f"  Total in tr_rexprt:      {n_rexprt:,}")
    print(f"  Predicted pathogenic:    {n_path:,}")
    con.close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--db",    default="data/trs.db")
    args = ap.parse_args()
    ingest(args.input, args.db)