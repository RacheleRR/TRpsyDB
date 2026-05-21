"""
Step 3 — parse TRExplorer TSV → populate trs + tr_established + tr_clinvar_functional
Run: python 03_ingest_trexplorer.py --tsv data/raw/TR_catalog.*.tsv --db data/trs.db --stripy data/raw/stripy_loci.json
"""
import sqlite3, gzip, argparse, csv, json
from pathlib import Path
from collections import defaultdict

SOURCES_INCLUDE = {
    "TRExplorerV1:KnownDiseaseAssociatedLoci",
    "TRExplorerV2:KnownDiseaseAssociatedLociV2",
    "TRExplorerV2:KnownFunctionalVNTRs",
    "TRExplorerV2:ClinvarIndelsThatAreTRs2025",
    "TRExplorerV1:Illumina174kPolymorphicTRs",
    "TRExplorerV1:PolymorphicTRsInT2TAssemblies",
    "TRExplorerV2:PolymorphicTRsInT2TAssembliesV2",
    "TRExplorerV2:AdottoTRsFromDanzi2025",
    "TRExplorerV2:VamosV3",
    "TRExplorerV2:Manigbas2024",
    "TRExplorerV2:Tanudisastro2025",
    "TRExplorerV2:Sulovari2021",
    "TRExplorerV2:Mukamel2021",
    "TRExplorerV2:Annear2021",
    "TRExplorerV2:Garg2021",
}

# Only these → tr_established (Mendelian disease with STRipy ranges)
SOURCES_MENDELIAN = {
    "TRExplorerV1:KnownDiseaseAssociatedLoci",
    "TRExplorerV2:KnownDiseaseAssociatedLociV2",
}

# These → tr_clinvar_functional
SOURCES_FUNCTIONAL = {
    "TRExplorerV2:KnownFunctionalVNTRs",
    "TRExplorerV2:ClinvarIndelsThatAreTRs2025",
}

GARG_SOURCE = "TRExplorerV2:Garg2021"

def safe_str(v):
    if v is None: return None
    if isinstance(v, list): v = ",".join(str(x) for x in v)
    s = str(v).strip()
    return s if s and s != "." else None

def safe_float(v):
    try: return float(v) if v and str(v).strip() not in (".", "") else None
    except: return None

def safe_int(v):
    try: return int(float(v)) if v and str(v).strip() not in (".", "") else None
    except: return None

def build_stripy_index(stripy_path):
    with open(stripy_path) as f:
        raw = json.load(f)
    by_coord = {}
    by_gene  = {}
    for locus_id, entry in raw.items():
        coords = entry.get("LocationCoordinates", {})
        hg38   = coords.get("hg38")
        gene   = entry.get("Gene") or locus_id
        if hg38:
            by_coord[hg38] = entry
        by_gene[gene] = entry
    print(f"  STRipy: {len(by_coord)} by coord, {len(by_gene)} by gene")
    return by_coord, by_gene

def extract_stripy_fields(entry):
    """Returns list of dicts — one per disease (some loci have multiple)."""
    diseases = entry.get("Diseases") or {}
    gene     = entry.get("Gene")
    motif    = entry.get("Motif")
    results  = []
    for dis_code, dis in diseases.items():
        normal   = dis.get("NormalRange") or {}
        intermed = dis.get("IntermediateRange") or {}
        results.append({
            "disease_name":   safe_str(dis.get("DiseaseName")),
            "gene_name":      safe_str(gene),
            "inheritance":    safe_str(dis.get("Inheritance")),
            "normal_min":     safe_int(normal.get("Min")),
            "normal_max":     safe_int(normal.get("Max")),
            "premut_min":     safe_int(intermed.get("Min")),
            "premut_max":     safe_int(intermed.get("Max")),
            "pathogenic_min": safe_int(dis.get("PathogenicCutoff")),
            "repeat_unit":    safe_str(motif),
            "pmid":           None,
        })
    return results

def ingest(tsv_path, db_path, stripy_path):

    # 1. STRipy index
    print("Loading STRipy data ...")
    by_coord, by_gene = build_stripy_index(stripy_path)

    # 2. Read and group by LocusId
    print(f"\nReading {tsv_path} ...")
    loci  = defaultdict(list)
    total = kept = 0

    opener = gzip.open if str(tsv_path).endswith(".gz") else open
    with opener(tsv_path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            total += 1
            if row["Source"] not in SOURCES_INCLUDE:
                continue
            kept += 1
            loci[row["LocusId"]].append(row)

    print(f"  Total rows:             {total:,}")
    print(f"  Rows in kept sources:   {kept:,}")
    print(f"  Unique loci (LocusId):  {len(loci):,}")

    # 3. Build rows
    print("\nBuilding rows ...")
    trs_rows      = []
    estab_rows    = []
    cvf_rows      = []
    stripy_matched = 0

    for i, (locus_id, rows) in enumerate(loci.items()):
        tr_id   = f"TR_{i+1:07d}"
        sources = sorted({r["Source"] for r in rows})
        ref     = rows[0]

        chrom       = safe_str(ref.get("chrom"))
        start       = safe_int(ref.get("start_0based"))
        end         = safe_int(ref.get("end_1based"))
        canon       = safe_str(ref.get("CanonicalMotif")) or safe_str(ref.get("ReferenceMotif"))
        ref_motif   = safe_str(ref.get("ReferenceMotif"))
        motif_size  = safe_int(ref.get("MotifSize"))
        n_repeats   = safe_int(ref.get("NumRepeatsInReference"))
        purity      = safe_float(ref.get("ReferenceRepeatPurity"))
        gene_name   = safe_str(ref.get("GencodeGeneName"))
        gene_region = safe_str(ref.get("GencodeGeneRegion"))
        gene_id     = safe_str(ref.get("GencodeGeneId"))
        mappability = safe_float(ref.get("FlanksAndLocusMappability"))
        noncoding   = safe_str(ref.get("NonCodingAnnotations"))
        af_illumina = safe_str(ref.get("AlleleFrequenciesFromIllumina174k"))
        af_t2t      = safe_str(ref.get("AlleleFrequenciesFromT2TAssemblies"))
        aou_median  = safe_float(ref.get("AoU1027_Median"))
        aou_stdev   = safe_float(ref.get("AoU1027_Stdev"))
        aou_max     = safe_float(ref.get("AoU1027_MaxAllele"))

        trs_rows.append((
            tr_id, locus_id,
            chrom, start, end,
            canon, ref_motif, motif_size,
            n_repeats, purity,
            gene_name, gene_region, gene_id,
            mappability, noncoding,
            af_illumina, af_t2t,
            aou_median, aou_stdev, aou_max,
            ",".join(sources), len(sources),
            int(GARG_SOURCE in sources),
            None, None, None,
        ))

        # tr_established — Mendelian disease loci only
        is_mendelian   = any(s in SOURCES_MENDELIAN for s in sources)
        is_functional  = any(s in SOURCES_FUNCTIONAL for s in sources)

        if is_mendelian:
            coord_str    = f"{chrom}:{start}-{end}"
            stripy_entry = by_coord.get(coord_str) or by_gene.get(gene_name)
            if stripy_entry:
                stripy_matched += 1
                diseases = extract_stripy_fields(stripy_entry)
            else:
                diseases = [{
                    "disease_name": None, "gene_name": gene_name,
                    "inheritance": None,
                    "normal_min": None, "normal_max": None,
                    "premut_min": None, "premut_max": None,
                    "pathogenic_min": None,
                    "repeat_unit": canon, "pmid": None,
                }]
            for dis in diseases:
                estab_rows.append((
                    tr_id,
                    dis["disease_name"],
                    dis["gene_name"] or gene_name,
                    dis["inheritance"],
                    dis["normal_min"], dis["normal_max"],
                    dis["premut_min"], dis["premut_max"],
                    dis["pathogenic_min"],
                    dis["repeat_unit"] or canon,
                    1, dis["pmid"],
                ))

        # tr_clinvar_functional — functional VNTRs and ClinVar loci
        if is_functional and not is_mendelian:
            for src in sources:
                if src in SOURCES_FUNCTIONAL:
                    locus_type = (
                        "functional_VNTR" if "KnownFunctional" in src
                        else "ClinVar"
                    )
                    cvf_rows.append((
                        tr_id, locus_type, gene_name, None, None, 3
                    ))
                    break

    # 4. Insert in chunks
    print("Inserting into database ...")
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA foreign_keys = ON")
    con.execute("PRAGMA journal_mode = WAL")
    con.execute("PRAGMA synchronous = NORMAL")
    cur = con.cursor()

    chunk = 50000
    for s in range(0, len(trs_rows), chunk):
        cur.executemany("""
            INSERT OR IGNORE INTO trs (
                tr_id, locus_id, chr, start, end,
                canonical_motif, reference_motif, motif_size_bp,
                num_repeats_ref, purity,
                gene_name, gene_region, gene_id,
                mappability, noncoding_annotations,
                af_illumina174k, af_t2t,
                aou_median, aou_stdev, aou_max,
                source_catalogs, n_sources, has_eVNTR,
                hg19_chr, hg19_start, hg19_end
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, trs_rows[s:s+chunk])
        con.commit()
        print(f"  trs: {min(s+chunk, len(trs_rows)):,} / {len(trs_rows):,}", end="\r")
    print()

    cur.executemany("""
        INSERT OR IGNORE INTO tr_established (
            tr_id, disease_name, gene_name, inheritance,
            normal_min, normal_max, premut_min, premut_max,
            pathogenic_min, repeat_unit, evidence_tier, pmid
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
    """, estab_rows)

    cur.executemany("""
        INSERT OR IGNORE INTO tr_clinvar_functional (
            tr_id, locus_type, gene_name, clinvar_id, evidence_note, evidence_tier
        ) VALUES (?,?,?,?,?,?)
    """, cvf_rows)

    con.commit()

    # 5. Summary
    n_trs   = cur.execute("SELECT COUNT(*) FROM trs").fetchone()[0]
    n_estab = cur.execute("SELECT COUNT(*) FROM tr_established").fetchone()[0]
    n_cvf   = cur.execute("SELECT COUNT(*) FROM tr_clinvar_functional").fetchone()[0]
    n_multi = cur.execute("SELECT COUNT(*) FROM trs WHERE n_sources > 1").fetchone()[0]
    n_evntr = cur.execute("SELECT COUNT(*) FROM trs WHERE has_eVNTR = 1").fetchone()[0]
    n_null  = cur.execute(
        "SELECT COUNT(*) FROM tr_established WHERE disease_name IS NULL"
    ).fetchone()[0]

    print(f"\nDone.")
    print(f"  trs rows:                  {n_trs:,}")
    print(f"  tr_established (Mendelian): {n_estab}")
    print(f"  tr_clinvar_functional:      {n_cvf}")
    print(f"  STRipy matched:             {stripy_matched} / {n_estab}")
    print(f"  multi-source TRs:           {n_multi:,}")
    print(f"  has_eVNTR = 1:              {n_evntr:,}")
    if n_null:
        print(f"\n  WARNING: {n_null} Mendelian loci unmatched to STRipy")
        print("  sqlite3 data/trs.db \"SELECT t.locus_id, t.gene_name")
        print("    FROM trs t JOIN tr_established e ON t.tr_id=e.tr_id")
        print("    WHERE e.disease_name IS NULL\"")
    con.close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv",    required=True)
    ap.add_argument("--db",     default="data/trs.db")
    ap.add_argument("--stripy", default="data/raw/stripy_loci.json")
    args = ap.parse_args()
    ingest(args.tsv, args.db, args.stripy)