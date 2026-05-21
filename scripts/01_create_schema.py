"""
Step 1 — create all empty tables in trs.db
Run: python 01_create_schema.py --db data/trs.db
"""
import sqlite3, argparse
from pathlib import Path

def create_schema(db_path):
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA foreign_keys = ON")
    cur = con.cursor()
    cur.executescript("""

    -- ── trs: one row per unique TR locus ─────────────────────────────────────
    CREATE TABLE IF NOT EXISTS trs (
        tr_id                   TEXT PRIMARY KEY,
        locus_id                TEXT UNIQUE,
        chr                     TEXT NOT NULL,
        start                   INTEGER NOT NULL,
        end                     INTEGER NOT NULL,
        canonical_motif         TEXT NOT NULL,
        reference_motif         TEXT,
        motif_size_bp           INTEGER,
        num_repeats_ref         INTEGER,
        purity                  REAL,
        gene_name               TEXT,
        gene_region             TEXT,
        gene_id                 TEXT,
        mappability             REAL,
        noncoding_annotations   TEXT,
        af_illumina174k         TEXT,
        af_t2t                  TEXT,
        aou_median              REAL,
        aou_stdev               REAL,
        aou_max                 REAL,
        source_catalogs         TEXT,
        n_sources               INTEGER DEFAULT 1,
        has_eVNTR               INTEGER DEFAULT 0,
        hg19_chr                TEXT,
        hg19_start              INTEGER,
        hg19_end                INTEGER
    );
    CREATE INDEX IF NOT EXISTS idx_trs_locus_id ON trs(locus_id);
    CREATE INDEX IF NOT EXISTS idx_trs_coord    ON trs(chr, start, end);
    CREATE INDEX IF NOT EXISTS idx_trs_gene     ON trs(gene_name);
    CREATE INDEX IF NOT EXISTS idx_trs_motif    ON trs(canonical_motif);

    -- ── tr_established: classical Mendelian TR expansion diseases ─────────────
    -- Source: KnownDiseaseAssociatedLoci V1+V2 only
    -- STRipy provides disease name, inheritance, repeat ranges
    CREATE TABLE IF NOT EXISTS tr_established (
        id              INTEGER PRIMARY KEY AUTOINCREMENT,
        tr_id           TEXT NOT NULL REFERENCES trs(tr_id),
        disease_name    TEXT,
        gene_name       TEXT,
        inheritance     TEXT,
        normal_min      INTEGER,
        normal_max      INTEGER,
        premut_min      INTEGER,
        premut_max      INTEGER,
        pathogenic_min  INTEGER,
        repeat_unit     TEXT,
        evidence_tier   INTEGER DEFAULT 1,
        pmid            TEXT
    );
    CREATE INDEX IF NOT EXISTS idx_estab_tr_id ON tr_established(tr_id);

    -- ── tr_clinvar_functional: functional VNTRs and ClinVar TR indels ─────────
    -- Source: KnownFunctionalVNTRs + ClinvarIndelsThatAreTRs2025
    CREATE TABLE IF NOT EXISTS tr_clinvar_functional (
        id              INTEGER PRIMARY KEY AUTOINCREMENT,
        tr_id           TEXT NOT NULL REFERENCES trs(tr_id),
        locus_type      TEXT NOT NULL,
        gene_name       TEXT,
        clinvar_id      TEXT,
        evidence_note   TEXT,
        evidence_tier   INTEGER DEFAULT 3
    );
    CREATE INDEX IF NOT EXISTS idx_cvf_tr_id ON tr_clinvar_functional(tr_id);

    -- ── tr_rexprt: RExPRT pathogenicity scores ────────────────────────────────
    CREATE TABLE IF NOT EXISTS tr_rexprt (
        tr_id           TEXT PRIMARY KEY REFERENCES trs(tr_id),
        ensembleScore   REAL,
        ensembleBinary  INTEGER,
        ensembleMax     REAL,
        SVM             REAL,
        XGB             REAL,
        pLi             REAL,
        loeuf           REAL,
        gc_content      REAL,
        per_g           REAL,
        per_c           REAL,
        per_a           REAL,
        per_t           REAL,
        eSTR            INTEGER,
        opReg           INTEGER,
        promoter        INTEGER,
        UTR_5           INTEGER,
        UTR_3           INTEGER,
        HG19_ID         TEXT
    );
    CREATE INDEX IF NOT EXISTS idx_rexprt_score ON tr_rexprt(ensembleScore);

    -- ── tr_expression: eVNTR / sc-eTR hits ───────────────────────────────────
    CREATE TABLE IF NOT EXISTS tr_expression (
        id               INTEGER PRIMARY KEY AUTOINCREMENT,
        tr_id            TEXT NOT NULL REFERENCES trs(tr_id),
        source           TEXT NOT NULL,
        qtl_type         TEXT NOT NULL,
        cell_type        TEXT,
        gene_affected    TEXT,
        gene_id          TEXT,
        beta             REAL,
        se               REAL,
        p_value          REAL,
        p_fdr            REAL,
        candidate_causal INTEGER DEFAULT 0,
        PIP              REAL,
        coloc_mQTL       INTEGER DEFAULT 0
    );
    CREATE INDEX IF NOT EXISTS idx_expr_tr_id ON tr_expression(tr_id);
    CREATE INDEX IF NOT EXISTS idx_expr_gene  ON tr_expression(gene_affected);

    -- ── tr_brain_qtl: TRxQTL brain QTL hits ──────────────────────────────────
    CREATE TABLE IF NOT EXISTS tr_brain_qtl (
        id          INTEGER PRIMARY KEY AUTOINCREMENT,
        tr_id       TEXT NOT NULL REFERENCES trs(tr_id),
        cohort      TEXT NOT NULL,
        tissue      TEXT NOT NULL,
        qtl_type    TEXT NOT NULL,
        gene_id     TEXT,
        beta        REAL,
        se          REAL,
        p_value     REAL,
        q_value     REAL
    );
    CREATE INDEX IF NOT EXISTS idx_bqtl_tr_id ON tr_brain_qtl(tr_id);
    CREATE INDEX IF NOT EXISTS idx_bqtl_gene  ON tr_brain_qtl(gene_id);

    -- ── tr_regulatory ─────────────────────────────────────────────────────────
    -- Column names match exactly what 07_annotate_regulatory.py outputs
    CREATE TABLE IF NOT EXISTS tr_regulatory (
        tr_id                    TEXT PRIMARY KEY REFERENCES trs(tr_id),

        -- EPD promoters (07_annotate_regulatory.py: annotate_epd)
        in_promoter              INTEGER DEFAULT 0,
        promoter_id              TEXT,

        -- ENCODE cCREs (07_annotate_regulatory.py: annotate_ccre)
        ccre_id                  TEXT,
        ccre_class               TEXT,
        in_any_ccre              INTEGER DEFAULT 0,
        brain_ccre               INTEGER DEFAULT 0,

        -- SEdb super enhancers (07_annotate_regulatory.py: annotate_sedb_se)
        in_any_SE                INTEGER DEFAULT 0,
        in_brain_SE              INTEGER DEFAULT 0,
        n_SE_tissues             INTEGER DEFAULT 0,
        SE_gene_overlap          TEXT,
        SE_gene_closest          TEXT,
        SE_gene_ABC              TEXT,
        SE_cas_value             REAL,

        -- SEdb typical enhancers (07_annotate_regulatory.py: annotate_sedb_te)
        in_any_TE                INTEGER DEFAULT 0,
        in_brain_TE              INTEGER DEFAULT 0,
        n_TE_tissues             INTEGER DEFAULT 0,

        -- TSS distances (07_annotate_regulatory.py: compute_tss_sj via R)
        dist_tss_bp              INTEGER,
        nearest_tss_gene         TEXT,

        -- SJ distances (07_annotate_regulatory.py: compute_tss_sj via R)
        dist_splice_bp           INTEGER
    );

    -- ── tr_phenotype_assoc: PheWAS hits ───────────────────────────────────────
    CREATE TABLE IF NOT EXISTS tr_phenotype_assoc (
        id                   INTEGER PRIMARY KEY AUTOINCREMENT,
        tr_id                TEXT NOT NULL REFERENCES trs(tr_id),
        source               TEXT NOT NULL,
        trait                TEXT,
        phenotype_code       TEXT,
        beta                 REAL,
        p_value              REAL,
        n_samples            INTEGER,
        psychiatric_relevant INTEGER DEFAULT 0,
        replicated_AoU       INTEGER DEFAULT 0
    );
    CREATE INDEX IF NOT EXISTS idx_pheno_tr_id ON tr_phenotype_assoc(tr_id);

    -- ── tr_psychiatric: psychiatric-focused TR associations ───────────────────
    CREATE TABLE IF NOT EXISTS tr_psychiatric (
        id                  INTEGER PRIMARY KEY AUTOINCREMENT,
        tr_id               TEXT NOT NULL REFERENCES trs(tr_id),
        phenotype           TEXT NOT NULL,
        phenotype_category  TEXT,
        study_type          TEXT,
        evidence_tier       INTEGER DEFAULT 3,
        pmid                TEXT,
        sc_eTR_overlap      INTEGER DEFAULT 0,
        gwas_coloc          INTEGER DEFAULT 0
    );
    CREATE INDEX IF NOT EXISTS idx_psych_tr_id ON tr_psychiatric(tr_id);

    -- ── gene_psychiatric_assoc: gene-level psychiatric evidence ───────────────
    CREATE TABLE IF NOT EXISTS gene_psychiatric_assoc (
        id               INTEGER PRIMARY KEY AUTOINCREMENT,
        gene_name        TEXT NOT NULL,
        gene_id          TEXT,
        phenotype        TEXT NOT NULL,
        source           TEXT NOT NULL,
        association_type TEXT,
        p_value          REAL,
        odds_ratio       REAL,
        pmid             TEXT,
        note             TEXT
    );
    CREATE INDEX IF NOT EXISTS idx_gpa_gene ON gene_psychiatric_assoc(gene_name);

    -- ── db_meta ───────────────────────────────────────────────────────────────
    CREATE TABLE IF NOT EXISTS db_meta (
        key   TEXT PRIMARY KEY,
        value TEXT
    );
    INSERT OR REPLACE INTO db_meta VALUES ('schema_version', '1.0');
    INSERT OR REPLACE INTO db_meta VALUES ('created', datetime('now'));
    """)

    con.commit()
    tables = cur.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    ).fetchall()
    print(f"Created: {db_path}")
    print(f"Tables ({len(tables)}):")
    for (t,) in tables:
        n = cur.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        print(f"  {t}: {n:,} rows")
    con.close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="data/trs.db")
    args = ap.parse_args()
    create_schema(args.db)