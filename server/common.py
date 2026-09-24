"""Workspace, settings and error helpers for the local server (not used by the browser build).

Every analysis the user runs is saved as an item under CB_HOME/work/<id>/: the input files, request.json (the
arguments, with file contents replaced by the saved file names) and result.json (what the UI, the PDF report
and the assistant read). CB_HOME defaults to ~/Library/Application Support/CloneBench.
"""
from __future__ import annotations

import base64
import json
import os
import re
import time
import traceback
import uuid
from pathlib import Path

from .core import UserFacingError  # noqa: F401  (re-exported: the one user-facing error class)

ROOT = Path(__file__).resolve().parent.parent
HOME = Path(os.environ.get("CB_HOME", Path.home() / "Library" / "Application Support" / "CloneBench"))
WORK = HOME / "work"
WORK.mkdir(parents=True, exist_ok=True)
SAVED_KINDS = ("construct", "primers", "sanger", "protein", "seqtools", "design")


def log_exc(where: str) -> str:
    """Print the current exception's traceback to the server log and return a short 'Type: msg' string."""
    import sys
    et, ev, _ = sys.exc_info()
    print(f"[{where}] non-fatal error:\n" + traceback.format_exc(), flush=True)
    return f"{et.__name__ if et else 'Error'}: {str(ev)[:160]}"


def _safe_name(n: str) -> str:
    n = Path(str(n)).name
    n = re.sub(r"[^\w.\- ()]+", "_", n).strip() or "file"
    return n[:120]


def _strip_files(obj, folder: Path, saved: list):
    """Write every {"name", "b64"} / {"name", "text"} file in the arguments into folder/input and replace it by a
    reference, so request.json stays small and the inputs are kept as the user gave them."""
    if isinstance(obj, dict):
        if "name" in obj and ("b64" in obj or ("text" in obj and len(obj) == 2)) and isinstance(obj.get("name"), str):
            p = folder / "input" / _safe_name(obj["name"])
            p.parent.mkdir(parents=True, exist_ok=True)
            data = base64.b64decode(obj["b64"]) if obj.get("b64") is not None else str(obj.get("text") or "").encode()
            p.write_bytes(data)
            saved.append(p.name)
            return {"name": obj["name"], "saved": f"input/{p.name}", "bytes": len(data)}
        out = {}
        for k, v in obj.items():
            if k == "record" and isinstance(v, dict) and "seq" in v:
                out[k] = {"name": v.get("name"), "length": v.get("length"), "from_item": v.get("from_item")}
            else:
                out[k] = _strip_files(v, folder, saved)
        return out
    if isinstance(obj, list):
        return [_strip_files(v, folder, saved) for v in obj]
    return obj


def save_item(kind: str, args: dict, result: dict) -> str:
    iid = time.strftime("%y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:4]
    d = WORK / iid
    d.mkdir(parents=True)
    saved = []
    req = {"kind": kind, "name": result.get("name"), "time": time.time(), "args": _strip_files(args, d, saved), "inputs": saved}
    (d / "request.json").write_text(json.dumps(req))
    result = dict(result)
    result["item"] = iid
    (d / "result.json").write_text(json.dumps(result))
    return iid


TOOL_NAMES = {"construct": "Construct", "primers": "Primers", "sanger": "Sanger", "protein": "Protein",
              "seqtools": "SequenceTools", "design": "PrimerDesign"}


def report_filename(iid: str, res: dict, ext: str = "pdf") -> str:
    """date_toolname_report — e.g. 2026-09-23_CloneBench-Sanger_report.pdf. The date is the day the analysis ran
    (from its request.json), so a report downloaded again later still carries the date of its contents."""
    try:
        t = json.loads((item_dir(iid) / "request.json").read_text()).get("time")
    except Exception:  # noqa: BLE001
        t = None
    day = time.strftime("%Y-%m-%d", time.localtime(t or time.time()))
    tool = TOOL_NAMES.get(res.get("kind"), str(res.get("kind") or "Analysis").title())
    return f"{day}_CloneBench-{tool}_report.{ext}"


def item_dir(iid: str) -> Path:
    p = (WORK / str(iid)).resolve()
    if not str(p).startswith(str(WORK.resolve())) or not (p / "result.json").exists():
        raise UserFacingError(f"No saved analysis with id {iid}.")
    return p


def load_item(iid: str) -> dict:
    r = json.loads((item_dir(iid) / "result.json").read_text())
    r["item"] = iid
    return r


def list_items(limit: int = 60, kind: str | None = None) -> list:
    rows = []
    for d in sorted(WORK.iterdir(), key=lambda x: x.name, reverse=True) if WORK.exists() else []:
        if not (d / "request.json").exists() or not (d / "result.json").exists():
            continue
        try:
            req = json.loads((d / "request.json").read_text())
            if kind and req.get("kind") != kind:
                continue
            rows.append({"item": d.name, "kind": req.get("kind"), "name": req.get("name"), "t": req.get("time", 0) * 1000})
        except Exception:  # noqa: BLE001 - a half-written item must not break the list
            continue
        if len(rows) >= limit:
            break
    return rows


def summary_of(iid: str) -> dict:
    try:
        return json.loads((WORK / iid / "result.json").read_text()).get("summary") or {}
    except Exception:  # noqa: BLE001
        return {}


# ---------------------------------------------------------------- settings (API key, model)
CONFIG = Path.home() / ".clone-bench" / "config.json"


def load_settings() -> dict:
    """Settings saved from the UI. The API key falls back to RNAseq Bench's, so one key serves every Bench app."""
    cfg = {}
    try:
        cfg = json.loads(CONFIG.read_text())
    except Exception:  # noqa: BLE001
        pass
    if not cfg.get("api_key"):
        try:
            k = json.loads((Path.home() / ".rnaseq-bench" / "config.json").read_text()).get("api_key")
            if k:
                cfg["api_key"] = k
        except Exception:  # noqa: BLE001
            pass
    return cfg


def save_settings(d: dict):
    CONFIG.parent.mkdir(parents=True, exist_ok=True)
    cur = {}
    try:
        cur = json.loads(CONFIG.read_text())
    except Exception:  # noqa: BLE001
        pass
    cur.update(d)
    CONFIG.write_text(json.dumps(cur))
    try:
        CONFIG.chmod(0o600)
    except Exception:  # noqa: BLE001
        pass


def api_key() -> str:
    return load_settings().get("api_key") or os.environ.get("ANTHROPIC_API_KEY") or ""
