"""Write examples/demo_reference.fasta (and .gb, with the read's ORF annotated) from examples/3730.ab1.

    .venv/bin/python tests/make_demo_reference.py

The reference is the read's own Mott-trimmed, high-quality basecalls with two planted edits: one substitution in a
high-quality stretch and one extra base (so the read shows a 1-bp deletion). It exists only so the Sanger tool has
something honest to show; it is not a real clone.
"""
import base64
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from server.core_sanger import build_demo_reference, demo_record  # noqa: E402
from server.core_seq import export_api  # noqa: E402

ab1 = {"name": "3730.ab1", "b64": base64.b64encode((ROOT / "examples" / "3730.ab1").read_bytes()).decode()}
d = build_demo_reference(ab1)
rec = demo_record(d)
head = (f">demo_reference 3730.ab1 basecalls (Mott-trimmed) with two planted edits: substitution at {d['substitution']['ref_pos1']} "
        f"(reference {d['substitution']['ref']}, read {d['substitution']['read']}); extra base {d['extra_base']['base']} at "
        f"{d['extra_base']['ref_pos1']} (the read shows a 1-bp deletion)")
seq = d["seq"]
(ROOT / "examples" / "demo_reference.fasta").write_text(head + "\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n")
gb = export_api({"record": rec, "format": "genbank"})["text"]
(ROOT / "examples" / "demo_reference.gb").write_text(gb)
print(json.dumps({k: v for k, v in d.items() if k != "seq"}, indent=1))
