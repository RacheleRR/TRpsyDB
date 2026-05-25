"""
TRpsyDB — FastAPI application
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, JSONResponse
from pathlib import Path
import os

from backend.routers import loci, analysis, meta, expression

BASE_DIR     = Path(__file__).parent.resolve()
FRONTEND_DIR = BASE_DIR / "frontend"

app = FastAPI(
    title="TRpsyDB",
    description="Neuropsychiatric Tandem Repeat Database",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(loci.router,       prefix="/api/loci",       tags=["loci"])
app.include_router(analysis.router,   prefix="/api/analysis",   tags=["analysis"])
app.include_router(meta.router,       prefix="/api/meta",       tags=["meta"])
app.include_router(expression.router, prefix="/api/expression", tags=["expression"])

def _serve(filename: str):
    p = FRONTEND_DIR / filename
    if p.exists():
        return FileResponse(str(p))
    return JSONResponse({"error": f"{filename} not found"}, status_code=404)

# ── Routes ────────────────────────────────────────────────────
@app.get("/",       include_in_schema=False)
async def serve_browse():   return _serve("browse.html")

@app.get("/browse", include_in_schema=False)
async def serve_browse2():  return _serve("browse.html")

@app.get("/locus",  include_in_schema=False)
async def serve_locus():    return _serve("locus.html")

# serve any other frontend file
@app.get("/{filename}.html", include_in_schema=False)
async def serve_html(filename: str): return _serve(f"{filename}.html")

if FRONTEND_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR)), name="static")

@app.get("/health", include_in_schema=False)
async def health():
    return {"status": "ok", "version": "1.0.0"}

@app.on_event("startup")
async def startup():
    from backend.db.database import get_db_stats
    print(f"[TRpsyDB] Frontend: {FRONTEND_DIR} (exists={FRONTEND_DIR.exists()})")
    try:
        stats = get_db_stats()
        print(f"[TRpsyDB] DB connected — {stats['n_loci']:,} loci, turso={stats['using_turso']}")
    except Exception as e:
        print(f"[TRpsyDB] DB connection error: {e}")
    print("[TRpsyDB] Ready →", os.environ.get("RENDER_EXTERNAL_URL", "http://localhost:8000"))