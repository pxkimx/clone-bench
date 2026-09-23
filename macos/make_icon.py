#!/usr/bin/env python3
"""Draw the app icon from the same geometry as the logo in web/index.html, and build AppIcon.icns.

    .venv/bin/python macos/make_icon.py

The coordinates below are the SVG's 32-unit viewBox (a rounded panel, the plasmid backbone ring, four feature arcs
and one cut site), scaled onto Apple's 1024 canvas, so the icon and the in-app mark cannot drift apart.
"""
from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
RES = HERE / "Clone Bench.app" / "Contents" / "Resources"
SS = 4
CANVAS, INSET, BOX = 1024, 100, 824
PANEL, RING, CUT = "#0d1b22", "#6a7b84", "#e4ecef"
# arcs on the r = 9 ring around (16, 16): (start degrees, end degrees, colour), 0° = 3 o'clock, clockwise
ARCS = [(-90, -24.3, "#00e0cf"), (16.8, 69.4, "#ffb020"), (109.4, 158.6, "#ff5c8a"), (198.0, 234.9, "#4ea3ff")]


def main():
    N = CANVAS * SS
    im = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    u = BOX * SS / 32.0
    o = INSET * SS
    P = lambda x, y: (o + x * u, o + y * u)  # noqa: E731
    d.rounded_rectangle([P(0.5, 0.5), P(31.5, 31.5)], radius=7.5 * u, fill=PANEL)
    # PIL draws an arc's width inwards from its bounding box, so grow the box by half the stroke to centre it on r = 9
    d.ellipse([P(7 - 0.65, 7 - 0.65), P(25 + 0.65, 25 + 0.65)], outline=RING, width=int(1.3 * u))
    w = int(3 * u)
    for a0, a1, c in ARCS:
        d.arc([P(7 - 1.5, 7 - 1.5), P(25 + 1.5, 25 + 1.5)], start=a0, end=a1, fill=c, width=w)
        for a in (a0, a1):                                   # round caps, as stroke-linecap="round"
            x, y = 16 + 9 * math.cos(math.radians(a)), 16 + 9 * math.sin(math.radians(a))
            d.ellipse([P(x - 1.5, y - 1.5), P(x + 1.5, y + 1.5)], fill=c)
    d.line([P(16, 3.2), P(16, 6.4)], fill=CUT, width=int(1.4 * u))
    for y in (3.2, 6.4):
        d.ellipse([P(16 - 0.7, y - 0.7), P(16 + 0.7, y + 0.7)], fill=CUT)
    im = im.resize((CANVAS, CANVAS), Image.LANCZOS)
    RES.mkdir(parents=True, exist_ok=True)
    im.save(RES / "AppIcon.png")
    if shutil.which("iconutil"):
        with tempfile.TemporaryDirectory() as td:
            iset = Path(td) / "AppIcon.iconset"
            iset.mkdir()
            for s in (16, 32, 128, 256, 512):
                im.resize((s, s), Image.LANCZOS).save(iset / f"icon_{s}x{s}.png")
                im.resize((2 * s, 2 * s), Image.LANCZOS).save(iset / f"icon_{s}x{s}@2x.png")
            subprocess.run(["iconutil", "-c", "icns", str(iset), "-o", str(RES / "AppIcon.icns")], check=True)
    print("wrote", RES / "AppIcon.png")


if __name__ == "__main__":
    main()
