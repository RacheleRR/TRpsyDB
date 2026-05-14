#!/bin/bash
echo "🧬 TRpsyDB starting up..."
echo "📦 Seeding database..."
python scripts/seed_classical_loci.py
echo "🚀 Starting server..."
uvicorn main:app --host 0.0.0.0 --port $PORT
