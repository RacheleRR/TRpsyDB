"""
Step 2 — fetch all loci from STRipy API and save to JSON
Run: python scripts/02_fetch_stripy.py --out data/raw/stripy_loci.json

STRipy API docs: https://stripy.org/docs
Endpoint used:  GET https://api.stripy.org/locus/{locus_id}
Returns: coordinates (hg38), motif, normal/intermediate/pathogenic ranges,
         disease name, gene, inheritance mode
"""
import json, time, argparse, urllib.request, urllib.error
from pathlib import Path

# All locus IDs in STRipy database (gene symbols used as IDs by their API)
# Source: https://stripy.org/database — 61 loci as of 2024
STRIPY_LOCUS_IDS = [
    "AFF2", "AFF3", "ARX", "ATN1", "ATXN1", "ATXN10", "ATXN2", "ATXN3",
    "ATXN7", "ATXN8OS", "BEAN1", "C9ORF72", "CACNA1A", "CBL", "CNBP",
    "COMP", "CSTB", "DAB1", "DMD", "DMPK", "EIF4A3", "FGF14", "FMR1",
    "FXN", "GIPC1", "GLS", "HOXA13a", "HOXA13b", "HOXA13c", "HOXD13",
    "HTT", "JPH3", "LRP12", "MARCHF6", "NIPA1", "NOP56", "NOTCH2NLC",
    "PABPN1", "PHOX2B", "PPP2R2B", "PRDM12", "PRNP", "RAPGEF2", "RFC1",
    "RILPL1", "RUNX2", "SAMD12", "SOX3", "STARD7", "TBP", "TCF4", "THAP11",
    "TNRC6A", "VLDLR", "VWA1", "WDR7", "XYLT1", "YEATS2", "ZIC2", "ZIC3",
    "ZNF713",
]

def fetch_locus(locus_id: str) -> dict | None:
    url = f"https://api.stripy.org/locus/{locus_id}"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=10) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        print(f"  HTTP {e.code} for {locus_id}")
        return None
    except Exception as e:
        print(f"  Error for {locus_id}: {e}")
        return None

def fetch_all(out_path: str):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    results = {}
    for i, locus_id in enumerate(STRIPY_LOCUS_IDS):
        print(f"  [{i+1}/{len(STRIPY_LOCUS_IDS)}] {locus_id} ...", end=" ")
        data = fetch_locus(locus_id)
        if data:
            results[locus_id] = data
            print("ok")
        else:
            print("FAILED")
        time.sleep(0.3)   # be polite to the API

    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved {len(results)}/{len(STRIPY_LOCUS_IDS)} loci to {out_path}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/raw/stripy_loci.json")
    args = ap.parse_args()
    fetch_all(args.out)