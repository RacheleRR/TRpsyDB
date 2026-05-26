#!/bin/bash
echo "🧬 TRpsyDB starting up..."
echo "🚀 Starting server..."
uvicorn main:app --host 0.0.0.0 --port $PORT
