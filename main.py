"""
TRpsyDB — FastAPI application
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, JSONResponse
from pathlib import Path
import os

from backend.db.database import init_db
from backend.routers import loci, analysis, meta, expression

# ── Paths ─────────────────────────────────────────────────────
# Use absolute path so it works both locally and on Render
BASE_DIR     = Path(__file__).parent.resolve()
FRONTEND_DIR = BASE_DIR / "frontend"
INDEX_HTML   = FRONTEND_DIR / "index.html"

# ── App ───────────────────────────────────────────────────────
app = FastAPI(
    title="TRpsyDB",
    description="Neuropsychiatric Tandem Repeat Database",
    version="0.1.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

# ── CORS ──────────────────────────────────────────────────────
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── API Routers ───────────────────────────────────────────────
app.include_router(loci.router,       prefix="/api/loci",       tags=["loci"])
app.include_router(analysis.router,   prefix="/api/analysis",   tags=["analysis"])
app.include_router(meta.router,       prefix="/api/meta",       tags=["meta"])
app.include_router(expression.router, prefix="/api/expression", tags=["expression"])

# ── Serve frontend HTML ───────────────────────────────────────
# This MUST come before the static mount so / returns index.html
@app.get("/", include_in_schema=False)
async def serve_frontend():
    if INDEX_HTML.exists():
        return FileResponse(str(INDEX_HTML))
    return JSONResponse({"error": f"index.html not found at {INDEX_HTML}"}, status_code=404)

# ── Serve static files (CSS, JS, assets inside frontend/) ────
if FRONTEND_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR)), name="static")
else:
    print(f"[WARNING] Frontend directory not found at {FRONTEND_DIR}")

# ── Health check ──────────────────────────────────────────────
@app.get("/health", include_in_schema=False)
async def health():
    return {"status": "ok", "version": "0.1.0"}

# ── Startup ───────────────────────────────────────────────────
@app.on_event("startup")
async def startup():
    print(f"[TRpsyDB] Base dir:     {BASE_DIR}")
    print(f"[TRpsyDB] Frontend dir: {FRONTEND_DIR} (exists={FRONTEND_DIR.exists()})")
    print(f"[TRpsyDB] index.html:   {INDEX_HTML} (exists={INDEX_HTML.exists()})")
    init_db()
    print("[TRpsyDB] Ready →", os.environ.get("RENDER_EXTERNAL_URL", "http://localhost:8000"))
