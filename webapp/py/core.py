"""Clone Bench core: the one entry point the UI, the assistant and the browser build all call.

    api(name, args) -> dict

Everything under server/core*.py is pure computation: Biopython, numpy and the standard library only, no
FastAPI, no file-system state. Arguments and results are JSON-serialisable dicts. The local server wraps this
with a workspace (history, PDF, assistant); the browser build (webapp/) runs these same files inside Pyodide,
which ships Biopython 1.85 rather than 1.88 — the few API differences are handled by shims in this package
(see `compat_aa_percent` in core_protein.py) rather than by dropping features.
"""
from __future__ import annotations

import base64
import platform
import re

import Bio


class UserFacingError(Exception):
    """An error whose message is written for the user (bad input, wrong file type) — shown without a traceback."""


# ---------------------------------------------------------------- result builder
class Result:
    """The shared result format (the same shape RNAseq Bench's result.json has).

    tiles    — headline numbers
    flags    — findings, level ok | warn | info; text may contain <b>, <i>, <code>
    sections — each with items: figures drawn in the browser from `data` (type "fig", `draw` names the
               renderer), tables (type "table") and interactive widgets (type "widget")
    methods  — [heading, paragraph] pairs ready for a methods section
    """

    def __init__(self, kind: str, name: str):
        self.r = {"kind": kind, "name": name, "tiles": [], "flags": [], "sections": [], "methods": [],
                  "versions": versions(), "summary": {}}
        self._sec = None

    def tile(self, value, label):
        self.r["tiles"].append({"value": value, "label": label})

    def flag(self, level, text):
        assert level in ("ok", "warn", "info")
        self.r["flags"].append({"level": level, "text": text})

    def section(self, sid, kicker, title, lede=""):
        self._sec = {"id": sid, "kicker": kicker, "title": title, "lede": lede, "items": []}
        self.r["sections"].append(self._sec)

    def fig(self, fid, title, draw, data, how="", yours="", wide=False, sub=""):
        self._sec["items"].append({"type": "fig", "id": fid, "title": title, "draw": draw, "data": data,
                                   "how": how, "yours": yours, "wide": wide, "sub": sub})

    def table(self, tid, title, columns, rows, note="", how="", yours="", wide=True):
        self._sec["items"].append({"type": "table", "id": tid, "title": title, "columns": columns, "rows": rows,
                                   "note": note, "how": how, "yours": yours, "wide": wide})

    def widget(self, wid, widget, title, data, how="", yours="", wide=True, sub=""):
        self._sec["items"].append({"type": "widget", "id": wid, "widget": widget, "title": title, "data": data,
                                   "how": how, "yours": yours, "wide": wide, "sub": sub})

    def method(self, head, text):
        self.r["methods"].append([head, text])

    def done(self):
        self.r["sections"] = [s for s in self.r["sections"] if s["items"]]
        return self.r


def versions() -> dict:
    v = {"biopython": Bio.__version__, "python": platform.python_version()}
    try:
        import numpy
        v["numpy"] = numpy.__version__
    except Exception:  # noqa: BLE001
        pass
    return v


# ---------------------------------------------------------------- small helpers
def file_bytes(f: dict) -> bytes:
    """A file sent by the UI: {"name": ..., "b64": ...} or {"name": ..., "text": ...}."""
    if not isinstance(f, dict):
        raise UserFacingError("No file was received.")
    if f.get("b64") is not None:
        try:
            return base64.b64decode(f["b64"])
        except Exception:  # noqa: BLE001
            raise UserFacingError(f"{f.get('name', 'The file')} could not be decoded.")
    if f.get("text") is not None:
        return f["text"].encode()
    raise UserFacingError("No file was received.")


def commas(v) -> str:
    return f"{int(round(float(v))):,}"


def esc(s) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def pl(n, word, plural=None) -> str:
    return f"{n:,} {word if n == 1 else (plural or word + 's')}"


def clean_name(s: str) -> str:
    s = re.sub(r"\\+", " ", str(s)).strip()
    return s[:80]


# ---------------------------------------------------------------- dispatcher
def api(name: str, args: dict | None = None) -> dict:
    """Run one named computation. Raises UserFacingError for bad input; anything else is a bug."""
    args = args or {}
    import importlib
    table = {
        # whole-tool analyses (the local app saves these as workspace items)
        "construct": ("core_seq", "construct"),
        "primers": ("core_primers", "primers"),
        "sanger": ("core_sanger", "sanger"),
        "protein": ("core_protein", "protein"),
        "seqtools": ("core_seq", "seqtools"),
        "design": ("core_design", "design"),
        # interactive follow-ups on a loaded construct
        "digest": ("core_seq", "digest_api"),
        "enzymes": ("core_seq", "enzymes_api"),
        "suggest_digest": ("core_seq", "suggest_api"),
        "orfs": ("core_seq", "orfs_api"),
        "translate_feature": ("core_seq", "translate_feature_api"),
        "export": ("core_seq", "export_api"),
        "read_reference": ("core_seq", "read_reference_api"),
        "demo_reference": ("core_sanger", "demo_reference_api"),
        "codon_tables": ("core_seq", "codon_tables_api"),
    }
    if name not in table:
        raise UserFacingError(f"Unknown request '{name}'.")
    mod, fn = table[name]
    return getattr(importlib.import_module(f".{mod}", __package__), fn)(args)
