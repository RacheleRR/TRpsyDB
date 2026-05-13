# TRpsyDB — Setup & Deployment Guide

## 1. Local setup (your machine)

```bash
# Clone after pushing to GitHub
git clone https://github.com/YOUR_USERNAME/TRpsyDB.git
cd TRpsyDB

# Install dependencies
pip install -r requirements.txt

# Run locally
uvicorn main:app --reload --port 8000

# Open http://localhost:8000
```

---

## 2. Push to GitHub (first time)

```bash
# Inside the TRpsyDB folder
git init
git add .
git commit -m "initial commit: TRpsyDB v0.1"

# Create repo on GitHub first (github.com → New repository → TRpsyDB)
# Then:
git remote add origin https://github.com/YOUR_USERNAME/TRpsyDB.git
git branch -M main
git push -u origin main
```

---

## 3. Deploy to Render (free, auto-deploys)

1. Go to https://render.com → sign up with GitHub
2. New → Web Service → Connect your TRpsyDB repo
3. Settings:
   - **Build command:** `pip install -r requirements.txt`
   - **Start command:** `uvicorn main:app --host 0.0.0.0 --port $PORT`
   - **Instance type:** Free
4. Click Deploy
5. Your live URL: `https://trpsydb.onrender.com`

Every `git push` to main → Render redeploys automatically.

---

## 4. Moving between machines (Ubuntu ↔ MacBook)

```bash
# On machine A: push changes
git add .
git commit -m "update"
git push

# On machine B: pull changes
git pull

# Database moves separately (not on GitHub)
# Copy psytrdb.sqlite manually or use the seed script:
python scripts/seed_db.py
```

---

## 5. Domain setup (optional, ~€12/year)

1. Buy `trpsydb.org` on Namecheap or Google Domains
2. In Render → Custom Domain → add `trpsydb.org`
3. In your domain registrar → add CNAME record pointing to Render URL

---

## 6. Adding your pipeline data

```bash
# Put your pipeline output in data/raw/
cp /path/to/ehdn_DBSCAN_annotated.tsv data/raw/
cp /path/to/rexprt_output.tsv data/raw/

# Run ingestion scripts (we build these next)
python scripts/ingest_pipeline.py
python scripts/ingest_gene_sets.py
python scripts/seed_classical_loci.py
```
