-- ============================================================
-- psyTRdb — Neuropsychiatric Tandem Repeat Database
-- Schema v0.1
-- ============================================================

PRAGMA foreign_keys = ON;

-- ------------------------------------------------------------
-- CORE: One row per TR locus
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS loci (
    id                  INTEGER PRIMARY KEY AUTOINCREMENT,

    -- Genomic coordinates (hg38)
    chrom               TEXT    NOT NULL,           -- e.g. chr3
    start               INTEGER NOT NULL,           -- 0-based
    end                 INTEGER NOT NULL,           -- 0-based
    locus_id            TEXT    UNIQUE NOT NULL,    -- e.g. chr3:63912684-63912714
    
    -- Repeat properties
    motif               TEXT    NOT NULL,           -- e.g. AAGG
    motif_length        INTEGER,                    -- derived: len(motif)
    reference_copies    REAL,                       -- copies in hg38
    normal_max          INTEGER,                    -- upper limit of normal range
    premutation_max     INTEGER,                    -- upper limit of premutation (if applicable)
    pathogenic_min      INTEGER,                    -- lower limit of pathogenic range (NULL if unknown)

    -- Gene context
    gene_symbol         TEXT,                       -- nearest gene (HGNC)
    ensembl_gene_id     TEXT,                       -- e.g. ENSG00000123456
    gene_region         TEXT,                       -- 5UTR / exon / intron / 3UTR / intergenic
    tss_distance        INTEGER,                    -- distance to nearest TSS (bp)
    splice_distance     INTEGER,                    -- distance to nearest splice junction (bp)
    strand              TEXT,                       -- + or -

    -- Regulatory / functional annotations
    rexprt_score        REAL,                       -- RExPRT pathogenicity score (0-1)
    encode_ccre         TEXT,                       -- ENCODE cCRE overlap (PLS/pELS/dELS/CTCF/DNase/NULL)
    roadmap_brain_state TEXT,                       -- Roadmap Epigenomics chromatin state in brain
    eqtl_evidence       TEXT,                       -- eQTL gene if known, else NULL

    -- Population genetics
    gnomad_mean_copies  REAL,                       -- gnomAD mean allele size
    gnomad_sd_copies    REAL,                       -- gnomAD SD
    gnomad_url          TEXT,                       -- direct link to gnomAD entry

    -- Psychiatric gene set membership (boolean flags — fast to query)
    in_schema           INTEGER DEFAULT 0,          -- 1 if gene in SCHEMA
    schema_or           REAL,                       -- SCHEMA odds ratio if available
    in_bipex            INTEGER DEFAULT 0,          -- 1 if gene in BipEx
    bipex_or            REAL,
    in_asd_gwas         INTEGER DEFAULT 0,
    in_mdd_gwas         INTEGER DEFAULT 0,
    in_adhd_gwas        INTEGER DEFAULT 0,
    in_scz_gwas         INTEGER DEFAULT 0,
    in_bpd_gwas         INTEGER DEFAULT 0,

    -- Gene constraint
    pli                 REAL,                       -- gnomAD pLI
    loeuf               REAL,                       -- gnomAD LOEUF

    -- Evidence summary (computed/updated when studies are added)
    evidence_tier       INTEGER DEFAULT 4,          -- 1=established 2=strong 3=moderate 4=exploratory 5=negative
    n_studies_tested    INTEGER DEFAULT 0,
    n_studies_positive  INTEGER DEFAULT 0,
    replication_status  TEXT DEFAULT 'untested',    -- untested/single/replicated/failed

    -- Associated phenotypes (comma-separated for simple cases)
    phenotypes          TEXT,                       -- e.g. "SCZ,BPD"

    -- Display
    common_name         TEXT,                       -- e.g. "RFC1 CANVAS repeat"
    notes               TEXT,                       -- free text, editorial notes

    -- Metadata
    created_at          TEXT DEFAULT (datetime('now')),
    updated_at          TEXT DEFAULT (datetime('now'))
);

-- Index for the most common queries
CREATE INDEX IF NOT EXISTS idx_loci_gene     ON loci(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_loci_chrom    ON loci(chrom, start, end);
CREATE INDEX IF NOT EXISTS idx_loci_tier     ON loci(evidence_tier);
CREATE INDEX IF NOT EXISTS idx_loci_schema   ON loci(in_schema);


-- ------------------------------------------------------------
-- STUDIES: Papers that have investigated TR loci
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS studies (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid            TEXT    UNIQUE,                 -- PubMed ID
    doi             TEXT,
    title           TEXT    NOT NULL,
    first_author    TEXT,
    year            INTEGER,
    journal         TEXT,
    cohort_size_cases    INTEGER,
    cohort_size_controls INTEGER,
    ancestry        TEXT,                           -- e.g. "European, East Asian"
    sequencing_type TEXT,                           -- WGS/WES/targeted/long-read
    tr_calling_tool TEXT,                           -- ExpansionHunter/EHdn/GangSTR/TRGT/etc
    phenotype       TEXT,                           -- SCZ/BPD/ASD/MDD/ADHD/mixed
    notes           TEXT,
    created_at      TEXT DEFAULT (datetime('now'))
);


-- ------------------------------------------------------------
-- EVIDENCE: Many-to-many between loci and studies
-- One row = one locus tested in one study
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS locus_evidence (
    id                  INTEGER PRIMARY KEY AUTOINCREMENT,
    locus_id            TEXT    NOT NULL REFERENCES loci(locus_id),
    study_id            INTEGER NOT NULL REFERENCES studies(id),

    -- What was found
    result              TEXT    NOT NULL,           -- 'positive' / 'negative' / 'trend'
    p_value             REAL,
    fdr                 REAL,
    odds_ratio          REAL,
    or_ci_lower         REAL,
    or_ci_upper         REAL,
    effect_direction    TEXT,                       -- 'expansion' / 'contraction' / 'both'
    
    -- Allele frequency in this study
    case_mean_copies    REAL,
    control_mean_copies REAL,
    case_af_outlier     REAL,                       -- freq of outlier alleles in cases
    control_af_outlier  REAL,

    -- Context
    analysis_type       TEXT,                       -- burden/single-locus/regression/etc
    covariates_used     TEXT,                       -- e.g. "PC1-10, sex, age"
    notes               TEXT,

    created_at          TEXT DEFAULT (datetime('now')),
    UNIQUE(locus_id, study_id)
);

CREATE INDEX IF NOT EXISTS idx_evidence_locus ON locus_evidence(locus_id);
CREATE INDEX IF NOT EXISTS idx_evidence_study ON locus_evidence(study_id);


-- ------------------------------------------------------------
-- GENE SETS: SCHEMA, BipEx, GWAS loci, SynGO, etc.
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS gene_sets (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    name        TEXT    UNIQUE NOT NULL,            -- e.g. "SCHEMA_FDR05"
    description TEXT,
    source_url  TEXT,
    version     TEXT,
    phenotype   TEXT,                               -- SCZ/BPD/ASD/etc
    n_genes     INTEGER,
    created_at  TEXT DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS gene_set_members (
    gene_set_id INTEGER NOT NULL REFERENCES gene_sets(id),
    gene_symbol TEXT    NOT NULL,
    ensembl_id  TEXT,
    odds_ratio  REAL,                               -- if available (SCHEMA, BipEx)
    p_value     REAL,
    notes       TEXT,
    PRIMARY KEY (gene_set_id, gene_symbol)
);

CREATE INDEX IF NOT EXISTS idx_gsm_gene ON gene_set_members(gene_symbol);


-- ------------------------------------------------------------
-- PATHWAY TERMS: Pre-computed pathway associations per gene
-- (populated from GO / Reactome / SynGO)
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS gene_pathways (
    gene_symbol TEXT    NOT NULL,
    pathway_id  TEXT    NOT NULL,                   -- e.g. GO:0007268
    pathway_name TEXT   NOT NULL,
    source      TEXT    NOT NULL,                   -- GO_BP / Reactome / SynGO
    PRIMARY KEY (gene_symbol, pathway_id)
);

CREATE INDEX IF NOT EXISTS idx_pathways_gene ON gene_pathways(gene_symbol);


-- ------------------------------------------------------------
-- USER ANALYSIS JOBS: Mode 2 — submitted TR lists
-- Ephemeral, cleared after 48h
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS analysis_jobs (
    id              TEXT    PRIMARY KEY,            -- UUID
    status          TEXT    DEFAULT 'pending',      -- pending/running/done/error
    submitted_at    TEXT    DEFAULT (datetime('now')),
    completed_at    TEXT,
    input_n_trs     INTEGER,
    result_json     TEXT                            -- full JSON result blob
);


-- ------------------------------------------------------------
-- METADATA: Key-value store for database versioning / stats
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS db_meta (
    key     TEXT PRIMARY KEY,
    value   TEXT
);

INSERT OR REPLACE INTO db_meta VALUES ('schema_version', '0.1');
INSERT OR REPLACE INTO db_meta VALUES ('created_at', datetime('now'));
INSERT OR REPLACE INTO db_meta VALUES ('n_loci', '0');
INSERT OR REPLACE INTO db_meta VALUES ('n_studies', '0');
