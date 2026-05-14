#!/usr/bin/env python3
"""
ingest_trexplorer.py
--------------------
Downloads and ingests the TRExplorer v2.0 catalog into TRpsyDB.

The catalog is a JSON-gz file (~5.6M loci). We ingest only the relevant
subset (polymorphic, non-homopolymer, mappability > 0.5) to keep the DB
manageable. You can adjust filters below.

Run from the TRpsyDB root:
    python scripts/ingest_trexplorer.py

Options:
    --download          Download the catalog first (needs ~1GB free space)
    --catalog PATH      Use a local catalog file instead of downloading
    --limit N           Only ingest first N loci (for testing)
    --min-mappability   Minimum mappability score (default: 0.5)
    --max-motif-size    Max motif size to include (default: 50)
    --no-homopolymers   Skip homopolymer (1bp motif) loci
    --reset             Clear existing TRExplorer loci before ingesting

Usage examples:
    # Download and ingest with defaults
    python scripts/ingest_trexplorer.py --download

    # Test run with first 10,000 loci
    python scripts/ingest_trexplorer.py --download --limit 10000

    # Use already-downloaded file
    python scripts/ingest_trexplorer.py --catalog data/raw/trexplorer_v2.json.gz
"""

import sqlite3
import json
import gzip
import argparse
import os
import sys
import urllib.request
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parent.parent
DB_PATH     = ROOT / "backend" / "db" / "psytrdb.sqlite"
DEFAULT_OUT = ROOT / "data" / "raw" / "trexplorer_v2.json.gz"

# ── Download URL ──────────────────────────────────────────────────────────────
# The EH-annotated JSON is the richest format (has RExPRT, gene, mappability...)
CATALOG_URL = (
    "https://github.com/broadinstitute/trexplorer-catalog/releases/download/v2.0/"
    "repeat_catalog_v2.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz"
)

# ── DB helpers ────────────────────────────────────────────────────────────────

def get_conn():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA cache_size=-65536")  # 64MB cache
    return conn


def ensure_tables(conn):
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS trs (
            tr_id           TEXT PRIMARY KEY,
            chrom           TEXT NOT NULL,
            start           INTEGER NOT NULL,
            end             INTEGER NOT NULL,
            motif           TEXT NOT NULL,
            motif_size      INTEGER,
            gene            TEXT,
            gene_region     TEXT,
            source_catalog  TEXT DEFAULT 'trexplorer_v2',
            rexprt_score    REAL,
            mappability     REAL,
            ref_purity      REAL,
            tss_dist        INTEGER,
            sj_dist         INTEGER,
            pli             REAL,
            loeuf           REAL,
            encode_ccre     TEXT,
            chromatin_state TEXT,
            is_established  INTEGER DEFAULT 0,
            created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );

        CREATE INDEX IF NOT EXISTS idx_trs_chrom_start ON trs(chrom, start);
        CREATE INDEX IF NOT EXISTS idx_trs_gene       ON trs(gene);
        CREATE INDEX IF NOT EXISTS idx_trs_rexprt     ON trs(rexprt_score);
    """)
    conn.commit()


# ── Parsing helpers ───────────────────────────────────────────────────────────

def parse_ref_region(ref_region):
    """Parse 'chr4:3074876-3074933' → (chr4, 3074876, 3074933)"""
    try:
        if isinstance(ref_region, list):
            ref_region = ref_region[0]
        chrom, coords = ref_region.split(":")
        start, end = coords.split("-")
        return chrom, int(start), int(end)
    except Exception:
        return None, None, None


def extract_motif(locus_structure):
    """Extract motif from '(CAG)*' → 'CAG'"""
    try:
        return locus_structure.strip("()*+").split(")*")[0].replace("(", "")
    except Exception:
        return locus_structure


def to_float(val):
    try:
        return float(val) if val is not None else None
    except (TypeError, ValueError):
        return None


# ── Main ingestion ────────────────────────────────────────────────────────────

def download_catalog(url: str, dest: Path):
    print(f"Downloading TRExplorer v2 catalog...")
    print(f"  URL: {url}")
    print(f"  Destination: {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)

    def progress(count, block_size, total_size):
        pct = count * block_size * 100 / total_size if total_size > 0 else 0
        mb  = count * block_size / 1024 / 1024
        print(f"\r  {mb:.0f} MB ({pct:.0f}%)  ", end="", flush=True)

    urllib.request.urlretrieve(url, dest, reporthook=progress)
    print(f"\n  Downloaded: {dest.stat().st_size / 1024 / 1024:.0f} MB")


def ingest_catalog(
    catalog_path: Path,
    limit: int = None,
    min_mappability: float = 0.5,
    max_motif_size: int = 50,
    no_homopolymers: bool = True,
    reset: bool = False,
):
    conn = get_conn()
    ensure_tables(conn)
    cur = conn.cursor()

    if reset:
        print("Resetting TRExplorer loci...")
        cur.execute("DELETE FROM trs WHERE source_catalog = 'trexplorer_v2'")
        conn.commit()

    print(f"\nIngesting from: {catalog_path}")
    print(f"  Filters: mappability≥{min_mappability}, motif_size≤{max_motif_size}, "
          f"{'no homopolymers' if no_homopolymers else 'including homopolymers'}")
    if limit:
        print(f"  Limit: {limit} loci (test mode)")

    inserted  = 0
    filtered  = 0
    errors    = 0
    batch     = []
    BATCH_SIZE = 5000

    opener = gzip.open if str(catalog_path).endswith(".gz") else open

    with opener(catalog_path, "rt", encoding="utf-8") as f:
        # The file is a JSON array — stream it line by line
        # Each locus is on its own line (one JSON object per line in the EH format)
        # but the file may also be a single large JSON array.
        # We handle both.

        raw = f.read(1)  # peek at first character
        f.seek(0)

        if raw == "[":
            # It's a JSON array — parse whole thing (memory intensive but correct)
            print("  Detected JSON array format — loading...")
            data = json.load(f)
            loci = data if isinstance(data, list) else data.get("loci", [])
        else:
            # Newline-delimited JSON
            loci = (json.loads(line) for line in f if line.strip())

        for locus in loci:
            if limit and (inserted + filtered) >= limit:
                break

            try:
                # ── Parse coordinates ─────────────────────────────────────────
                ref_region = locus.get("ReferenceRegion", "")
                chrom, start, end = parse_ref_region(ref_region)
                if chrom is None:
                    filtered += 1
                    continue

                # ── Motif ─────────────────────────────────────────────────────
                motif = extract_motif(locus.get("LocusStructure", ""))
                if not motif:
                    filtered += 1
                    continue

                motif_size = len(motif)

                # ── Filters ───────────────────────────────────────────────────
                if no_homopolymers and motif_size == 1:
                    filtered += 1
                    continue
                if motif_size > max_motif_size:
                    filtered += 1
                    continue

                mappability = to_float(locus.get("Mappability"))
                if mappability is not None and mappability < min_mappability:
                    filtered += 1
                    continue

                # ── Extract annotations ───────────────────────────────────────
                tr_id       = locus.get("LocusId") or f"{chrom}-{start}-{end}-{motif}"
                gene        = locus.get("Gene") or locus.get("GeneName")
                gene_region = locus.get("GeneRegion")
                rexprt      = to_float(locus.get("RExPRTScore") or locus.get("RExPRT"))
                ref_purity  = to_float(locus.get("RefPurity") or locus.get("BaseRepeatPurity"))
                source      = locus.get("Source", "trexplorer_v2")

                batch.append((
                    tr_id, chrom, start, end, motif, motif_size,
                    gene, gene_region, source,
                    rexprt, mappability, ref_purity
                ))

                inserted += 1

                # ── Flush batch ───────────────────────────────────────────────
                if len(batch) >= BATCH_SIZE:
                    cur.executemany("""
                        INSERT OR IGNORE INTO trs
                            (tr_id, chrom, start, end, motif, motif_size,
                             gene, gene_region, source_catalog,
                             rexprt_score, mappability, ref_purity)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """, batch)
                    conn.commit()
                    batch.clear()
                    print(f"\r  Inserted: {inserted:,}  Filtered: {filtered:,}  "
                          f"Errors: {errors}", end="", flush=True)

            except Exception as e:
                errors += 1
                if errors <= 5:
                    print(f"\n  Warning: {e} — locus: {locus.get('LocusId', '?')}")

    # Final flush
    if batch:
        cur.executemany("""
            INSERT OR IGNORE INTO trs
                (tr_id, chrom, start, end, motif, motif_size,
                 gene, gene_region, source_catalog,
                 rexprt_score, mappability, ref_purity)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, batch)
        conn.commit()

    conn.close()

    print(f"\n\n✅ Done.")
    print(f"   Inserted:  {inserted:,}")
    print(f"   Filtered:  {filtered:,}")
    print(f"   Errors:    {errors}")
    print(f"\n   DB: {DB_PATH}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ingest TRExplorer v2 catalog into TRpsyDB")
    parser.add_argument("--download",         action="store_true",
                        help="Download catalog from GitHub before ingesting")
    parser.add_argument("--catalog",          type=Path, default=DEFAULT_OUT,
                        help=f"Path to catalog JSON.gz (default: {DEFAULT_OUT})")
    parser.add_argument("--limit",            type=int, default=None,
                        help="Ingest only first N loci (useful for testing)")
    parser.add_argument("--min-mappability",  type=float, default=0.5,
                        help="Minimum mappability score (default: 0.5)")
    parser.add_argument("--max-motif-size",   type=int, default=50,
                        help="Max motif size in bp (default: 50)")
    parser.add_argument("--no-homopolymers",  action="store_true", default=True,
                        help="Skip 1-bp motif (homopolymer) loci (default: True)")
    parser.add_argument("--reset",            action="store_true",
                        help="Clear existing TRExplorer loci before ingesting")
    args = parser.parse_args()

    if args.download:
        download_catalog(CATALOG_URL, args.catalog)
    elif not args.catalog.exists():
        print(f"ERROR: Catalog not found at {args.catalog}")
        print("Run with --download to fetch it first.")
        sys.exit(1)

    ingest_catalog(
        catalog_path    = args.catalog,
        limit           = args.limit,
        min_mappability = args.min_mappability,
        max_motif_size  = args.max_motif_size,
        no_homopolymers = args.no_homopolymers,
        reset           = args.reset,
    )
