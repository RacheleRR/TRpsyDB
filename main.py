"""
psyTRdb — FastAPI application
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pathlib import Path

from backend.db.database import init_db, get_db_stats
from backend.routers import loci, analysis, meta, expression

# ── App ──────────────────────────────────────────────────────
app = FastAPI(
    title="TRpsyDB",
    description="Neuropsychiatric Tandem Repeat Database",
    version="0.1.0",
)

# ── CORS (allow frontend to call API during local dev) ───────
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],        # tighten this before production
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Routers ──────────────────────────────────────────────────
app.include_router(loci.router,     prefix="/api/loci",     tags=["loci"])
app.include_router(analysis.router, prefix="/api/analysis", tags=["analysis"])
app.include_router(meta.router,       prefix="/api/meta",       tags=["meta"])
app.include_router(expression.router, prefix="/api/expression", tags=["expression"])

# ── Serve frontend ───────────────────────────────────────────
FRONTEND_DIR = Path(__file__).parent / "frontend"
if FRONTEND_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR)), name="static")

@app.get("/", include_in_schema=False)
def serve_frontend():
    return FileResponse(str(FRONTEND_DIR / "index.html"))

# ── Startup ──────────────────────────────────────────────────
@app.on_event("startup")
def startup():
    init_db()
    print("[psyTRdb] Ready")
