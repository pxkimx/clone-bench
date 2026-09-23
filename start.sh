#!/usr/bin/env bash
# Clone Bench — first run creates a Python environment (~1 min), later runs start in seconds.
set -e
cd "$(dirname "$0")"
PY=$(command -v python3.12 || command -v python3.13 || command -v python3.11 || command -v python3.10 || command -v python3)
if [ ! -d .venv ]; then
  echo "Creating Python environment with $PY…"
  "$PY" -m venv .venv
  .venv/bin/pip install --upgrade pip -q
  .venv/bin/pip install -r requirements.txt
fi
PORT=${PORT:-8768}
( sleep 3; open "http://localhost:$PORT" 2>/dev/null || xdg-open "http://localhost:$PORT" 2>/dev/null ) &
exec .venv/bin/python -m uvicorn server.app:app --port "$PORT"
