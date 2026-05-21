"""
Step 4 — export ALL trs loci → RExPRT input TSV
Run: python 04_prepare_rexprt_input.py --db data/trs.db --out data/processed/rexprt_input_all.tsv
Then: bash rexptr_standalone.sh data/processed/rexprt_input_all.tsv
"""
import sqlite3, argparse, csv
from pathlib import Path

def export(db_path, out_path):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    rows = cur.execute("""
        SELECT tr_id, chr, start, end, canonical_motif, gene_name
        FROM trs
        ORDER BY chr, start
    """).fetchall()
    con.close()

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        # contig/start/end/motif = what rexptr_standalone.sh expects
        # tr_id = our key to match back after RExPRT
        w.writerow(["contig", "start", "end", "motif", "gene", "tr_id"])
        for tr_id, chrom, start, end, motif, gene in rows:
            w.writerow([chrom, start, end, motif or "", gene or "", tr_id])

    print(f"Written {len(rows):,} loci → {out_path}")
    print(f"\nNext:")
    print(f"  cd ~/test_TREX/TREX")
    print(f"  bash rexptr_standalone.sh {Path(out_path).resolve()}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--db",  default="data/trs.db")
    ap.add_argument("--out", default="data/processed/rexprt_input_all.tsv")
    args = ap.parse_args()
    export(args.db, args.out)