"""In-app Claude assistant. Streams a tool-using conversation over SSE ({type: text|tool|tool_result|done|error}),
the same protocol as RNAseq Bench and MassSpec Bench. Every tool calls the same core functions the UI uses."""
from __future__ import annotations

import json
import os
import threading
from typing import Iterator

from . import core
from .common import UserFacingError, list_items, load_item, load_settings

DEFAULT_MODEL = "claude-sonnet-5"
FALLBACK_MODELS = ["claude-sonnet-5", "claude-opus-5-5", "claude-haiku-4-5-20251001"]


def available_models() -> dict:
    cfg = load_settings()
    key = cfg.get("api_key") or os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        return {"ok": False, "error": "no key", "models": FALLBACK_MODELS}
    try:
        import anthropic
        client = anthropic.Anthropic(api_key=key)
        ids = [m.id for m in client.models.list(limit=100)]
        ids = [i for i in ids if i.startswith("claude")]
        ids = sorted(ids, key=lambda i: (0 if "sonnet" in i else 1 if "opus" in i else 2))
        return {"ok": True, "models": ids or FALLBACK_MODELS}
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "error": f"{type(e).__name__}: {str(e)[:200]}", "models": FALLBACK_MODELS}


def _resolve_model(client, wanted: str) -> str:
    """Use the configured model if the key can see it; otherwise the newest Sonnet (models get retired)."""
    try:
        ids = [m.id for m in client.models.list(limit=100)]
    except Exception:  # noqa: BLE001
        return wanted
    if wanted in ids:
        return wanted
    for pat in ("sonnet", "opus", "haiku"):
        for i in ids:
            if pat in i:
                return i
    return wanted


SYSTEM = """You are the cloning assistant built into Clone Bench, a local app for bench cloning work built on Biopython.
It has five tools the user sees: Construct (plasmid map, features and CDS checks, restriction analysis with commercially
available enzymes, virtual digest and agarose gel, ORFs), Primers (Tm by nearest-neighbour thermodynamics under a stated
buffer, 3' clamp, repeats, a complementarity heuristic for dimers and hairpins, virtual PCR on a template), Sanger
(trace viewer, Mott trimming, alignment to a reference, every difference called with a confidence, mixed peaks, which
parts of the construct are verified), Protein (ProtParam properties, A280 to concentration, hydropathy, tags and protease
sites) and Sequence tools.

You talk to a bench biologist looking at a result page. Explain in plain language what a result means for their clone and
what to do next. Quote numbers from tool results; never invent a cut position, a fragment size, a Tm or a mutation — call
a tool. Coordinates are 1-based; restriction cut positions are the first base after the cut on the top strand. Be honest
about limits: the dimer/hairpin check is a complementarity heuristic without free energies (primer3 or vendor tools give
ΔG); a Sanger difference at Phred < 20 may not be real; the virtual gel is an illustration; a digest compared with "the
construct minus a feature" stands in for the empty vector, which may differ. Vendor Tm calculators (NEB, Thermo) use
their own parameters and will differ by a few °C. Be concise: short paragraphs or bullet lists, the answer first."""

ITEM = {"type": "string", "description": "workspace item id (from list_items); omit to use the item the user is viewing"}
TOOLS = [
    {"name": "list_items", "description": "Recent analyses in the workspace: id, kind (construct, primers, sanger, protein, seqtools) and name.",
     "input_schema": {"type": "object", "properties": {"kind": {"type": "string"}}}},
    {"name": "get_item_result", "description": "The result of one analysis: tiles, findings, each figure's plain-language reading, table rows, and for a construct its feature list. Read this before interpreting anything.",
     "input_schema": {"type": "object", "properties": {"item": ITEM}}},
    {"name": "digest", "description": "Virtual digest of a construct (or any item with a sequence). enzymes: a list of enzyme names for one lane, or a list of lists for several lanes. Optional agarose (0.7, 1 or 2 %) and compare_feature (a feature name: also digests the construct with that feature deleted, a stand-in for the empty vector). Returns fragments, co-migrating bands and a verdict.",
     "input_schema": {"type": "object", "properties": {"item": ITEM, "enzymes": {"type": "array", "items": {}}, "agarose": {"type": "number"},
                                                       "compare_feature": {"type": "string"}}, "required": ["enzymes"]}},
    {"name": "find_enzymes", "description": "Restriction enzymes (commercially available) filtered: 'single', 'double', 'all', 'none' (does not cut), 'outside' (cuts once, not inside feature), 'inside' (cuts only inside feature), 'flank' (single cutters nearest each side of feature, with distances). feature = feature name for the last three.",
     "input_schema": {"type": "object", "properties": {"item": ITEM, "filter": {"type": "string"}, "feature": {"type": "string"}}, "required": ["filter"]}},
    {"name": "suggest_digest", "description": "Search one- or two-enzyme digests whose fragments are all resolvable and well spread; with a feature, only digests that cut inside it and differ from the construct without it (the empty-vector stand-in).",
     "input_schema": {"type": "object", "properties": {"item": ITEM, "feature": {"type": "string"}, "agarose": {"type": "number"}}}},
    {"name": "primer_check", "description": "Analyse primers: Tm (nearest-neighbour, GC, Wallace), GC, 3' clamp, repeats, dimer/hairpin heuristic, pair analysis, and a virtual PCR when template_item is given. primers: text, one per line, 'name sequence'. buffer: {preset: 'taq'|'hifi'|'custom', Na, K, Tris, Mg, dNTPs (mM), primer_nM, DMSO (%)}.",
     "input_schema": {"type": "object", "properties": {"primers": {"type": "string"}, "buffer": {"type": "object"}, "template_item": {"type": "string"}}, "required": ["primers"]}},
    {"name": "sanger_summary", "description": "Summary of a Sanger analysis item: per read coverage, identity, every difference with position, quality, confidence and consequence, mixed peaks, and which features are verified.",
     "input_schema": {"type": "object", "properties": {"item": ITEM}}},
    {"name": "translate_feature", "description": "DNA and protein sequence of a feature (by name) in an item's sequence, with the CDS check for CDS features.",
     "input_schema": {"type": "object", "properties": {"item": ITEM, "feature_name": {"type": "string"}}, "required": ["feature_name"]}},
    {"name": "protein_properties", "description": "ProtParam properties of a protein sequence (or DNA, translated): MW, pI, charge at 7.4, extinction coefficient, Abs 0.1 %, GRAVY, instability, tags and protease sites.",
     "input_schema": {"type": "object", "properties": {"seq": {"type": "string"}}, "required": ["seq"]}},
    {"name": "find_orfs", "description": "Six-frame ORF scan of an item's sequence (circular sequences across the origin). min_aa default 75, table = NCBI genetic code (default 11).",
     "input_schema": {"type": "object", "properties": {"item": ITEM, "min_aa": {"type": "integer"}, "table": {"type": "integer"}}}},
]


def _record(item: str | None, current: str | None) -> tuple[dict, dict]:
    iid = item or current
    if not iid:
        raise UserFacingError("No item given and the user is not viewing one. Call list_items.")
    r = load_item(iid)
    rj = r.get("record")
    if not rj:
        raise UserFacingError(f"Item {iid} ({r.get('kind')}) has no DNA sequence attached.")
    return r, rj


def _slim(r: dict) -> dict:
    out = {k: r.get(k) for k in ("item", "kind", "name", "tiles", "flags", "methods", "summary")}
    secs = []
    for s in r.get("sections", []):
        items = []
        for it in s["items"]:
            x = {"type": it["type"], "id": it.get("id"), "title": it.get("title"), "yours": it.get("yours")}
            if it["type"] == "table":
                x["columns"] = it["columns"]
                x["rows"] = it["rows"][:25]
                x["n_rows"] = len(it["rows"])
            if it["type"] == "fig" and it.get("draw") == "mono":
                x["text"] = it["data"].get("text", "")[:1500]
            if it["type"] == "widget" and it.get("widget") == "digest":
                x["suggested_digests"] = it["data"].get("suggestions")
            items.append(x)
        secs.append({"id": s["id"], "title": s["title"], "items": items})
    out["sections"] = secs
    rj = r.get("record")
    if rj:
        out["sequence"] = {"name": rj["name"], "length": rj["length"], "circular": rj["circular"],
                           "features": [{"name": f["name"], "type": f["type"], "start": f["start0"] + 1, "end": f["end0"],
                                         "strand": f["strand"], "wraps_origin": f["wraps"]} for f in rj.get("features", [])]}
    return out


def run_tool(name: str, inp: dict, current: str | None) -> str:
    try:
        if name == "list_items":
            return json.dumps(list_items(30, inp.get("kind")))
        if name == "get_item_result":
            return json.dumps(_slim(load_item(inp.get("item") or current)))[:60000]
        if name == "digest":
            _, rj = _record(inp.get("item"), current)
            enz = inp.get("enzymes") or []
            lanes = enz if enz and all(isinstance(e, list) for e in enz) else [enz]
            d = core.api("digest", {"record": rj, "lanes": lanes, "agarose": inp.get("agarose") or 1.0,
                                    "compare_feature": inp.get("compare_feature")})
            for lane in d["lanes"] + d["compare"]:
                for b in lane["bands"]:
                    b.pop("y", None)
            return json.dumps(d)[:40000]
        if name == "find_enzymes":
            _, rj = _record(inp.get("item"), current)
            d = core.api("enzymes", {"record": rj, "filter": inp.get("filter") or "single", "feature": inp.get("feature")})
            for k in ("rows", "upstream", "downstream"):
                if k in d:
                    d[k] = [{x: r[x] for x in ("enzyme", "cuts", "positions", "site", "overhang", "features", "distance") if x in r}
                            for r in d[k][:80]]
            if "names" in d:
                d["names"] = d["names"][:400]
            return json.dumps(d)[:40000]
        if name == "suggest_digest":
            _, rj = _record(inp.get("item"), current)
            return json.dumps(core.api("suggest_digest", {"record": rj, "feature": inp.get("feature"), "agarose": inp.get("agarose") or 1.0}))
        if name == "primer_check":
            b = inp.get("buffer") or {}
            args = {"text": inp["primers"], "preset": b.get("preset") or "taq", "buffer": {k: v for k, v in b.items() if k != "preset"}}
            if inp.get("template_item"):
                _, rj = _record(inp["template_item"], None)
                args["template"] = rj
            r = core.api("primers", args)
            r.pop("record", None)
            return json.dumps(_slim(r) | {"buffer": r.get("buffer"), "pcr_products": (r.get("pcr") or {}).get("products")})[:50000]
        if name == "sanger_summary":
            r = load_item(inp.get("item") or current)
            if r.get("kind") != "sanger":
                return f"error: item {r.get('item')} is a {r.get('kind')} analysis, not Sanger."
            s = _slim(r)
            s.pop("methods", None)
            return json.dumps(s)[:60000]
        if name == "translate_feature":
            _, rj = _record(inp.get("item"), current)
            return json.dumps(core.api("translate_feature", {"record": rj, "feature": inp["feature_name"]}))
        if name == "protein_properties":
            r = core.api("protein", {"text": inp["seq"]})
            s = _slim(r)
            s.pop("methods", None)
            return json.dumps(s)[:30000]
        if name == "find_orfs":
            _, rj = _record(inp.get("item"), current)
            d = core.api("orfs", {"record": rj, "min_aa": inp.get("min_aa") or 75, "table": inp.get("table") or 11})
            d["orfs"] = [{k: o[k] for k in ("frame", "strand", "start0", "end0", "aa", "match", "wraps")} | {"start": o["start0"] + 1}
                         for o in d["orfs"][:60]]
            return json.dumps(d)
        return f"unknown tool {name}"
    except UserFacingError as e:
        return f"error: {e}"
    except Exception as e:  # noqa: BLE001
        from .common import log_exc
        log_exc(f"agent tool {name}")
        return f"error: {type(e).__name__}: {e}"


CONV: dict[str, list] = {}
LOCK = threading.Lock()


def chat(conv_id: str, user_text: str, item: str | None) -> Iterator[str]:
    """Yields SSE 'data:' lines: {type: text|tool|tool_result|done|error}."""
    try:
        import anthropic
    except ImportError:
        yield _sse({"type": "error", "text": "The anthropic package is missing: pip install anthropic"})
        return
    cfg = load_settings()
    key = cfg.get("api_key") or os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        yield _sse({"type": "error", "text": "No API key. Open Settings (gear icon) and paste an Anthropic API key."})
        return
    client = anthropic.Anthropic(api_key=key)
    model = _resolve_model(client, cfg.get("model") or DEFAULT_MODEL)
    with LOCK:
        hist = CONV.setdefault(conv_id, [])
    ctx = ""
    if item:
        try:
            r = load_item(item)
            ctx = (f"\n\nThe user is viewing item '{item}' ({r['kind']}): {r['name']}. Findings:\n"
                   + "\n".join(f"- [{f['level']}] {f['text']}" for f in r["flags"][:20]))
        except Exception:  # noqa: BLE001
            ctx = ""
    hist.append({"role": "user", "content": user_text})
    for _ in range(12):
        try:
            with client.messages.stream(model=model, max_tokens=4000, system=SYSTEM + ctx, tools=TOOLS, messages=hist) as stream:
                for ev in stream:
                    if ev.type == "content_block_delta" and getattr(ev.delta, "type", "") == "text_delta":
                        yield _sse({"type": "text", "text": ev.delta.text})
                msg = stream.get_final_message()
        except Exception as e:  # noqa: BLE001
            m = str(e)
            if "authentication" in m.lower() or "invalid x-api-key" in m.lower():
                m = "The API key was rejected. Open Settings & API key and paste a valid key from console.anthropic.com."
            elif "credit" in m.lower() or "billing" in m.lower():
                m = "Anthropic says this key has no credit. Add a prepaid balance at console.anthropic.com → Billing."
            yield _sse({"type": "error", "text": f"{type(e).__name__}: {m}"})
            hist.pop()
            return
        hist.append({"role": "assistant", "content": [b.model_dump() for b in msg.content]})
        tool_uses = [b for b in msg.content if b.type == "tool_use"]
        if not tool_uses:
            break
        results = []
        for tu in tool_uses:
            yield _sse({"type": "tool", "name": tu.name, "input": tu.input})
            out = run_tool(tu.name, tu.input, item)
            yield _sse({"type": "tool_result", "name": tu.name, "text": out[:1500]})
            results.append({"type": "tool_result", "tool_use_id": tu.id, "content": out})
        hist.append({"role": "user", "content": results})
    yield _sse({"type": "done"})


def _sse(d: dict) -> str:
    return f"data: {json.dumps(d)}\n\n"


def reset(conv_id: str):
    with LOCK:
        CONV.pop(conv_id, None)
