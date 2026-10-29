"""The page's router must know every tool the sidebar links to.

"Design primers" was added to the sidebar in 0.2.0 but not to the list in route(), so `#tool=design` — a reload, the
Back button, a link — showed the home page instead of the tool. Clicking worked, which is why it went unnoticed.
"""
import re
from pathlib import Path

HTML = (Path(__file__).resolve().parent.parent / "web" / "index.html").read_text(encoding="utf-8")


def test_every_sidebar_tool_is_routable():
    linked = set(re.findall(r'href="#tool=(\w+)"', HTML))
    assert linked, "no tool links found in the sidebar"
    m = re.search(r"function route\(\)\{.*?\[([^\]]+)\]\.includes\(t\[1\]\)", HTML, re.S)
    assert m, "route() no longer has a list of tools — update this test"
    routable = set(re.findall(r'"(\w+)"', m.group(1)))
    assert linked <= routable, f"route() does not handle: {sorted(linked - routable)}"
    for tool in linked:
        assert f'id="v_{tool}"' in HTML, f"no view for tool {tool}"
