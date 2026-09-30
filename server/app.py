"""Clone Bench web server.  Run:  uvicorn server.app:app --port 8768

The UI calls one function, api(name, args). Here that is POST /api/run/<name> with the arguments as JSON; the
server runs the same core.api the browser build runs in Pyodide, and — for the five whole-tool analyses — saves
the result as a workspace item (history, PDF report, assistant).
"""
from __future__ import annotations

import asyncio
import os
import platform
import threading

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import agent, core, lifecycle
from .common import (HOME, ROOT, SAVED_KINDS, UserFacingError, item_dir, list_items, load_item, load_settings, log_exc,
                     report_filename, save_item, save_settings, summary_of)

VERSION = (ROOT / "VERSION").read_text().strip() if (ROOT / "VERSION").exists() else "dev"
app = FastAPI(title="Clone Bench", version=VERSION)
lifecycle.start()  # the server stops when the last page is closed, or when the launcher is gone
print(f"Clone Bench {VERSION} — workspace in {HOME} — python {platform.python_version()} — "
      f"biopython {core.versions()['biopython']}", flush=True)


@app.post("/api/run/{name}")
async def run(name: str, request: Request):
    try:
        args = await request.json()
    except Exception:  # noqa: BLE001
        raise HTTPException(400, "The request was not valid JSON.")
    if not isinstance(args, dict):
        raise HTTPException(400, "The request must be a JSON object.")
    try:
        # the computation is CPU-bound and short; run it off the event loop
        import anyio
        res = await anyio.to_thread.run_sync(core.api, name, args)
    except UserFacingError as e:
        raise HTTPException(400, str(e))
    except HTTPException:
        raise
    except Exception as e:  # noqa: BLE001
        msg = log_exc(f"api {name}")
        raise HTTPException(500, f"Something went wrong inside Clone Bench ({msg}). The details are in the server log.") from e
    if name in SAVED_KINDS and not args.get("_nosave"):
        try:
            res["item"] = save_item(name, args, res)
        except Exception:  # noqa: BLE001 - failing to save must not lose the result on screen
            log_exc("save item")
    return JSONResponse(res)


@app.get("/api/items")
def items(limit: int = 60, kind: str | None = None):
    rows = list_items(limit, kind)
    for r in rows:
        r["summary"] = summary_of(r["item"])
    return rows


@app.get("/api/items/{iid}")
def item(iid: str):
    try:
        return JSONResponse(load_item(iid))
    except UserFacingError as e:
        raise HTTPException(404, str(e))


@app.get("/api/items/{iid}/files/{path:path}")
def item_file(iid: str, path: str):
    try:
        d = item_dir(iid)
    except UserFacingError as e:
        raise HTTPException(404, str(e))
    p = (d / path).resolve()
    if not str(p).startswith(str(d)) or not p.is_file():
        raise HTTPException(404)
    return FileResponse(p, filename=p.name)


class ReportIn(BaseModel):
    figures: dict = {}
    widget_tables: dict = {}


@app.post("/api/items/{iid}/report")
def report(iid: str, body: ReportIn):
    """PDF of a saved analysis. The page sends PNG snapshots of its figures (they are drawn in the browser)."""
    from .report import build_pdf
    try:
        d = item_dir(iid)
        res = load_item(iid)
    except UserFacingError as e:
        raise HTTPException(404, str(e))
    out = d / "report.pdf"
    try:
        build_pdf(res, out, body.figures, body.widget_tables)
    except Exception as e:  # noqa: BLE001
        msg = log_exc("pdf")
        raise HTTPException(500, f"The PDF report could not be written ({msg}).") from e
    return FileResponse(out, media_type="application/pdf", filename=report_filename(iid, res))


# ---------------------------------------------------------------- settings, assistant
class Settings(BaseModel):
    api_key: str | None = None
    model: str | None = None


@app.get("/api/settings")
def get_settings():
    cfg = load_settings()
    key = cfg.get("api_key") or os.environ.get("ANTHROPIC_API_KEY") or ""
    return {"version": VERSION, "has_key": bool(key), "key_hint": (key[:7] + "…" + key[-4:]) if key else "",
            "model": cfg.get("model") or agent.DEFAULT_MODEL, "home": str(HOME), "biopython": core.versions()["biopython"]}


@app.post("/api/settings")
def set_settings(body: Settings):
    save_settings({k: v for k, v in body.model_dump().items() if v})
    return get_settings()


@app.get("/api/models")
def list_models():
    return agent.available_models()


@app.post("/api/quit")
def quit_server():
    """Stop the server (sidebar 'Quit'). Saved analyses stay on disk."""
    threading.Timer(0.5, lambda: lifecycle.quit_now("asked to quit")).start()
    return {"ok": True}


# Closing the window closes the app. An open page holds the connection below, which drops the moment
# the window goes away even if the browser closes the tab without telling the page; it also sends a
# heartbeat, and says so on the way out when it gets the chance. See server/lifecycle.py for the rule.
@app.get("/api/window")
async def window_channel(request: Request, page: str = ""):
    """Held open for as long as a page is open. Nothing is sent over it but keep-alives."""
    async def stream():
        token = lifecycle.window_opened(page)
        try:
            while not await request.is_disconnected():
                yield b": open\n\n"
                await asyncio.sleep(lifecycle.CHECK_S)
        finally:
            lifecycle.window_closed(token)

    return StreamingResponse(stream(), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})


@app.post("/api/alive")
def page_alive(page: str = ""):
    lifecycle.ping(page or "page")
    return {"ok": True}


@app.post("/api/closing")
def page_closing(page: str = ""):
    lifecycle.closing(page or "page")
    return {"ok": True}


class ChatIn(BaseModel):
    conv: str
    message: str
    item: str | None = None


@app.post("/api/agent/chat")
def agent_chat(body: ChatIn):
    return StreamingResponse(agent.chat(body.conv, body.message, body.item), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})


@app.post("/api/agent/reset")
def agent_reset(body: ChatIn):
    agent.reset(body.conv)
    return {"ok": True}


app.mount("/examples", StaticFiles(directory=ROOT / "examples"), name="examples")
app.mount("/", StaticFiles(directory=ROOT / "web", html=True), name="web")
