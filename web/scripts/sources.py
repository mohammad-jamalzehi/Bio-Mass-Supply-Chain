"""Publish browsable snapshots of cited Python and Markdown files."""
from html import escape
from pathlib import Path

WEB = Path(__file__).resolve().parents[1]
ROOT = WEB.parent
SOURCES = (
    "optimization/phase3.py", "optimization/annual.py",
    "validation/phase1.py", "tests/test_phase2.py", "tests/test_phase3.py",
    "docs/phase-4-annual.md",
)

for name in SOURCES:
    original = (ROOT / name).read_text(encoding="utf-8")
    rows = "\n".join(
        f'<span id="L{n}"><a href="#L{n}">{n:4}</a> {escape(line)}</span>'
        for n, line in enumerate(original.splitlines(), 1)
    )
    title = escape(name)
    destination = WEB / "public/research/sources" / (name + ".html")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        f"""<!doctype html><html lang="en"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="author" content="Mohamad Jamalzehi, Kosar Ghazali"><title>{title}</title>
<style>
body{{font:16px system-ui,sans-serif;margin:1.5rem;color:#193943;background:#f8f8f3}}
pre{{overflow:auto;padding:1rem;border:1px solid #dce4db;background:white;line-height:1.5}}
pre span{{display:block}}pre span:target{{background:#dbeee4}}
a{{color:#116d61}}
</style></head><body><a href="../../../">Research website</a><h1>{title}</h1>
<p>Project authors: Mohamad Jamalzehi and Kosar Ghazali</p>
<pre><code>{rows}</code></pre></body></html>
""", encoding="utf-8",
    )
