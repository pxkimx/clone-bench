"""PDF report of a saved analysis (reportlab), in the MassSpec Bench style.

Figures are drawn in the browser, so the page sends PNG snapshots of them (and the current rows of its interactive
tables) when you ask for the report; the text — tiles, findings, how/yours, tables, methods — comes from
result.json. A figure the page did not send is listed by title with its text, never silently dropped.
"""
from __future__ import annotations

import base64
import io
import re
from datetime import datetime
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (Image as RLImage, KeepTogether, PageBreak, Paragraph, Preformatted, SimpleDocTemplate,
                                Spacer, Table, TableStyle)

INK = colors.HexColor("#13212b")
MUTED = colors.HexColor("#5b6b75")
ACCENT = colors.HexColor("#0f766e")
WARN = colors.HexColor("#b45309")
LINE = colors.HexColor("#d8e0e4")

KIND_LABEL = {"construct": "Construct", "primers": "Primers", "sanger": "Sanger verification", "protein": "Protein",
              "seqtools": "Sequence tools"}

S = {
    "title": ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=22, leading=26, textColor=INK),
    "kicker": ParagraphStyle("kicker", fontName="Helvetica-Bold", fontSize=8, leading=10, textColor=ACCENT, spaceAfter=2),
    "h1": ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=16, leading=20, textColor=INK, spaceAfter=4),
    "h2": ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=11, leading=14, textColor=INK, spaceBefore=4, spaceAfter=2),
    "body": ParagraphStyle("body", fontName="Helvetica", fontSize=9, leading=12.5, textColor=INK, alignment=TA_LEFT),
    "muted": ParagraphStyle("muted", fontName="Helvetica", fontSize=8.5, leading=11.5, textColor=MUTED),
    "sub": ParagraphStyle("sub", fontName="Helvetica-Oblique", fontSize=7.5, leading=10, textColor=MUTED),
    "cell": ParagraphStyle("cell", fontName="Helvetica", fontSize=7, leading=9, textColor=INK),
    "mono": ParagraphStyle("mono", fontName="Courier", fontSize=6.5, leading=8, textColor=INK),
    "tile_v": ParagraphStyle("tv", fontName="Helvetica-Bold", fontSize=14, leading=17, textColor=INK),
    "tile_l": ParagraphStyle("tl", fontName="Helvetica", fontSize=7.5, leading=9, textColor=MUTED),
}
TRANS = {"−": "-", "→": "->", "±": "+/-", "–": "-", "—": "-", "×": "x", "≈": "~", "²": "2", "⁺": "+", "⁻": "-", "¹": "1",
         "′": "'", "↓": "|", "≥": "&gt;=", "≤": "&lt;=", "ε": "e", "Δ": "d", "·": "-", "…": "...", "“": '"', "”": '"',
         "’": "'", "₂": "2", "µ": "u", "β": "beta", "✓": "ok"}


def clean(html) -> str:
    html = str(html or "")
    # a bare "<" in prose ("Q < 20") is not a tag; left in, the tag-stripper below ate everything from it to the next
    # ">" — an opening <font> included — and ReportLab failed on the orphaned </font>, so no Sanger report could be written
    html = re.sub(r"<(?![a-zA-Z/])", "&lt;", html)
    html = re.sub(r"<code>(.*?)</code>", r"<font face='Courier'>\1</font>", html)
    html = re.sub(r"<(?!/?(b|i|font)\b)[^>]+>", "", html)
    for a, b in TRANS.items():
        html = html.replace(a, b)
    return html


def para(html, style):
    """A Paragraph from app HTML; text ReportLab still cannot parse goes in as plain text rather than losing the report."""
    try:
        return Paragraph(clean(html), style)
    except Exception:  # noqa: BLE001
        from xml.sax.saxutils import escape
        return Paragraph(escape(plain(re.sub(r"<[^>]+>", "", str(html or "")))), style)


def plain(s) -> str:
    s = str(s or "")
    for a, b in TRANS.items():
        s = s.replace(a, b.replace("&gt;", ">").replace("&lt;", "<"))
    return s


def _png(data_url: str):
    if not data_url or "," not in data_url:
        return None
    try:
        raw = base64.b64decode(data_url.split(",", 1)[1])
        img = ImageReader(io.BytesIO(raw))
        return raw, img.getSize()
    except Exception:  # noqa: BLE001
        return None


def _table(cols, rows, width, max_rows=60):
    data = [[Paragraph(f"<b>{clean(c)}</b>", S["cell"]) for c in cols]] + \
           [[para(v, S["cell"]) for v in r] for r in rows[:max_rows]]
    t = Table(data, repeatRows=1)
    t.setStyle(TableStyle([("LINEBELOW", (0, 0), (-1, 0), 0.8, INK), ("LINEBELOW", (0, 1), (-1, -1), 0.3, LINE),
                           ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
                           ("VALIGN", (0, 0), (-1, -1), "TOP")]))
    return t


def build_pdf(res: dict, out: Path, figures: dict | None = None, widget_tables: dict | None = None):
    figures = figures or {}
    widget_tables = widget_tables or {}
    W, H = A4
    margin = 16 * mm
    width = W - 2 * margin
    story = [Paragraph("CLONE BENCH · " + KIND_LABEL.get(res["kind"], res["kind"]).upper() + " REPORT", S["kicker"]),
             para(res["name"], S["title"]), Spacer(1, 3),
             Paragraph(datetime.now().strftime("Generated %d %B %Y, %H:%M"), S["muted"]), Spacer(1, 10)]
    cells = [[para(str(t["value"]), S["tile_v"]), para(t["label"], S["tile_l"])] for t in res["tiles"]]
    rows = [[Table([[c[0]], [c[1]]], style=[("LEFTPADDING", (0, 0), (-1, -1), 0)]) for c in cells[i:i + 3]] for i in range(0, len(cells), 3)]
    if rows:
        for r in rows:
            while len(r) < 3:
                r.append("")
        t = Table(rows, colWidths=[width / 3] * 3)
        t.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.6, LINE), ("INNERGRID", (0, 0), (-1, -1), 0.6, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 7), ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
                               ("LEFTPADDING", (0, 0), (-1, -1), 9), ("VALIGN", (0, 0), (-1, -1), "TOP")]))
        story += [t, Spacer(1, 12)]
    story.append(Paragraph("What your data says", S["h1"]))
    fl = []
    for f in res["flags"]:
        tag = {"ok": ("OK", ACCENT), "warn": ("CHECK", WARN), "info": ("NOTE", MUTED)}[f["level"]]
        fl.append([Paragraph(f"<font color='#{tag[1].hexval()[2:]}'><b>{tag[0]}</b></font>", S["cell"]),
                   para(f["text"], S["body"])])
    if fl:
        t = Table(fl, colWidths=[16 * mm, width - 16 * mm])
        t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, -2), 0.4, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                               ("LEFTPADDING", (0, 0), (-1, -1), 0)]))
        story.append(t)

    for sec in res["sections"]:
        story += [PageBreak(), para(str(sec["kicker"] or "").upper(), S["kicker"]), para(sec["title"], S["h1"]),
                  para(sec.get("lede"), S["muted"]), Spacer(1, 8)]
        for it in sec["items"]:
            block = [para(it.get("title"), S["h2"])]
            if it["type"] in ("fig", "widget"):
                png = _png(figures.get(it["id"]))
                if png:
                    raw, (w, h) = png
                    dw = width if it.get("wide") else width * 0.75
                    dw = min(dw, w * 0.5)
                    dh = dw * h / w
                    if dh > 170 * mm:
                        dh = 170 * mm
                        dw = dh * w / h
                    block += [Spacer(1, 3), RLImage(io.BytesIO(raw), width=dw, height=dh, hAlign="LEFT"), Spacer(1, 4)]
                elif it["type"] == "fig" and it.get("draw") == "mono":
                    block.append(Preformatted(plain(it["data"].get("text", ""))[:6000], S["mono"]))
                elif it["type"] == "fig" and it.get("draw") in ("alignment",):
                    pass
                else:
                    block.append(Paragraph("<i>(interactive figure — open the analysis in Clone Bench to see it)</i>", S["sub"]))
                for wt in widget_tables.get(it["id"], []):
                    block += [Spacer(1, 3), para(wt.get("title"), S["sub"]), _table(wt["columns"], wt["rows"], width)]
            elif it["type"] == "table":
                block.append(_table(it["columns"], it["rows"], width))
                if it.get("note"):
                    block.append(para(it["note"], S["sub"]))
            if it.get("how"):
                block.append(Paragraph("<b>How to read it.</b> " + clean(it["how"]), S["muted"]))
            if it.get("yours"):
                block.append(Paragraph("<b>In your data.</b> " + clean(it["yours"]), S["body"]))
            block.append(Spacer(1, 12))
            story.append(KeepTogether(block) if len(block) < 8 else block[0])
            if len(block) >= 8:
                story += block[1:]

    story += [PageBreak(), Paragraph("METHODS", S["kicker"]), Paragraph("How this analysis was computed", S["h1"])]
    for h, p in res["methods"]:
        story += [para(h, S["h2"]), para(p, S["body"])]
    if res.get("versions"):
        story += [Spacer(1, 8), Paragraph("Software versions", S["h2"]),
                  Paragraph(", ".join(f"{k} {v}" for k, v in res["versions"].items()), S["muted"])]
    story += [Spacer(1, 10), Paragraph("Interpretations are generated automatically from the numbers above and are meant to "
                                       "guide, not replace, expert review.", S["sub"])]

    def deco(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(margin, 10 * mm, "Clone Bench · created by Paul H. Kim, Ph.D.")
        canvas.drawRightString(W - margin, 10 * mm, f"{doc.page}")
        canvas.setStrokeColor(ACCENT)
        canvas.setLineWidth(2)
        canvas.line(margin, H - 10 * mm, margin + 18 * mm, H - 10 * mm)
        canvas.restoreState()

    doc = SimpleDocTemplate(str(out), pagesize=A4, leftMargin=margin, rightMargin=margin, topMargin=16 * mm, bottomMargin=16 * mm,
                            title=f"Clone Bench report — {plain(res['name'])}", author="Paul H. Kim, Ph.D.")
    doc.build(story, onFirstPage=deco, onLaterPages=deco)
    return out
