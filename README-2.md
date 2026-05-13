# TRpsyDB — Neuropsychiatric Tandem Repeat Database

[![License: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/Database-CC%20BY--NC%204.0-lightgrey.svg)](data/LICENSE_DATABASE.txt)
[![GTEx v10](https://img.shields.io/badge/Expression-GTEx%20v10-green.svg)](https://gtexportal.org)
[![bioRxiv](https://img.shields.io/badge/Preprint-coming%20soon-orange.svg)]()

## License

| Component | License |
|---|---|
| All code (FastAPI, scripts, pipeline) | [MIT](LICENSE) — free for any use |
| Curated database content | [CC BY-NC 4.0](data/LICENSE_DATABASE.txt) — free for academic use, commercial use requires permission |

Third-party data sources (GTEx, SCHEMA, BipEx, gnomAD) retain their own licenses — see [data/LICENSE_DATABASE.txt](data/LICENSE_DATABASE.txt).

A curated, evidence-tiered web database of tandem repeat loci with relevance 
to neuropsychiatric disorders (SCZ, BPD, ASD, MDD, ADHD).

## Project structure

```
psyTRdb/
├── main.py                    # FastAPI app entry point
├── requirements.txt
├── frontend/
│   └── index.html             # Single-file frontend (no framework needed)
├── backend/
│   ├── db/
│   │   ├── schema.sql         # Full database schema
│   │   ├── database.py        # DB connection + helpers
│   │   └── psytrdb.sqlite     # SQLite database (created on first run)
│   └── routers/
│       ├── loci.py            # Mode 1: gene/locus search
│       ├── analysis.py        # Mode 2: TR list enrichment
│       └── meta.py            # DB stats
├── data/
│   ├── raw/                   # Put your pipeline output files here
│   └── processed/             # Cleaned/formatted data ready for import
└── scripts/                   # Data ingestion scripts (next step)
```

## Local development

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Run the server
uvicorn main:app --reload --port 8000

# 3. Open http://localhost:8000 in your browser
```

## API endpoints

| Endpoint | Description |
|---|---|
| `GET /api/loci/search?q=RFC1` | Search by gene or locus ID |
| `GET /api/loci/{locus_id}` | Full locus report card |
| `GET /api/loci/gene/{gene}` | All TRs near a gene |
| `GET /api/loci/?tier=2&phenotype=SCZ` | Browse with filters |
| `POST /api/analysis/enrich` | Mode 2: TR list enrichment |
| `GET /api/meta/` | Database stats |

## Next steps

1. **Add your pipeline data** → put TSV files in `data/raw/`, 
   then we write ingestion scripts in `scripts/`
2. **Add classical loci** → HTT, FMR1, C9orf72, RFC1, ATXN... (Tier 1)
3. **Mine the literature** → extract loci from PsychENCODE-era TR papers
4. **Add TSS BED file** → for full Mode 2 proximity annotation
5. **Deploy to Render** → push to GitHub → connect to Render.com free tier

## Deployment (Render.com)

1. Push this repo to GitHub
2. Go to render.com → New → Web Service
3. Connect your GitHub repo
4. Build command: `pip install -r requirements.txt`
5. Start command: `uvicorn main:app --host 0.0.0.0 --port $PORT`
6. Done — free tier, auto-deploys on every git push
