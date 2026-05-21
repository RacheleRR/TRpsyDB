"""
Step 5 — clean RExPRT output → ready for ingestion
Run: python 05_clean_rexprt_output.py \
        --rexprt  ~/test_TREX/TREX/results/rexprt/rexprt_input_all_rex_input_TRsAnnotated_RExPRTscores.txt \
        --input   data/processed/rexprt_input_all.tsv \
        --out     data/processed/rexprt_cleaned.tsv

What this does (replaces your R script entirely):
  1. Reads raw RExPRT output
  2. Creates HG19_ID = chr#start#end
  3. Drops intermediate model columns
  4. Merges with original input TSV to recover tr_id and gene
  5. Writes clean TSV ready for step 6
"""
import csv, argparse
from pathlib import Path

# Columns to drop — same as your R script
COLS_DROP = {
    "chr", "start", "end", "id",
    "eSh0", "eSh1", "eSh2", "eSh3", "eSh4", "eSh5",
    "eTr0", "eTr1", "eTr2", "eTr3", "eTr4", "eTr5",
    "eH", "eW", "eS6", "eS", "J",
    "eX0", "eX1R", "eX2", "eX3", "eX4", "eX5",
    "tissue_simple_No_expression", "tissue_simple_Nervous_System",
    "tissue_simple_Other",
    "region_intergenic", "region_intron", "region_exon",
    "location_First", "location_Middle", "location_Last",
}

def clean(rexprt_path, input_path, out_path):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    # 1. Index our input TSV: RExPRT ID → tr_id
    # RExPRT ID format: motif#contig#start#end (built by rexptr_standalone.sh)
    print(f"Loading input TSV: {input_path}")
    id_to_tr = {}
    with open(input_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rex_id = f"{row['motif']}#{row['contig']}#{row['start']}#{row['end']}"
            id_to_tr[rex_id] = row.get("tr_id", "")
    print(f"  {len(id_to_tr):,} loci indexed")

    # 2. Read RExPRT output
    print(f"Loading RExPRT output: {rexprt_path}")
    with open(rexprt_path) as f:
        raw_rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"  {len(raw_rows):,} rows")
    if raw_rows:
        print(f"  Columns: {list(raw_rows[0].keys())}")

    # 3. Clean each row
    cleaned   = []
    unmatched = 0

    for row in raw_rows:
        # HG19_ID = chr#start#end (matches your R code)
        hg19_id = f"{row.get('chr','')}#{row.get('start','')}#{row.get('end','')}"

        # recover tr_id from ID column
        rex_id = row.get("ID", "")
        tr_id  = id_to_tr.get(rex_id, "")
        if not tr_id:
            unmatched += 1

        out_row = {"HG19_ID": hg19_id, "tr_id": tr_id}
        for col, val in row.items():
            if col not in COLS_DROP:
                out_row[col] = val
        cleaned.append(out_row)

    if not cleaned:
        print("ERROR: no rows to write")
        return

    # 4. Write
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cleaned[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(cleaned)

    print(f"\nDone.")
    print(f"  Cleaned rows: {len(cleaned):,}")
    print(f"  Unmatched:    {unmatched}")
    print(f"  Output:       {out_path}")
    print(f"\nNext:")
    print(f"  python 06_ingest_rexprt.py --input {out_path} --db data/trs.db")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rexprt", required=True,
                    help="Raw RExPRT output (*_TRsAnnotated_RExPRTscores.txt)")
    ap.add_argument("--input",  required=True,
                    help="Original input TSV from step 4")
    ap.add_argument("--out",    default="data/processed/rexprt_cleaned.tsv")
    args = ap.parse_args()
    clean(args.rexprt, args.input, args.out)