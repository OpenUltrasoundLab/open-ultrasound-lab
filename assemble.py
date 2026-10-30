#!/usr/bin/env python3
"""Assemble module fragments into the shared shell.
Usage: python3 assemble.py b1..b5     -> build/part1-foundations.html ... build/part5-modes-capstones.html
                                         (Part I..V, rail<N>/home<N>; default b1)
       python3 assemble.py m3         -> build/m3-only.html, the module with its part's rail/home
"""
import re, sys, pathlib
root = pathlib.Path(__file__).parent
BATCH = {1: list(range(1, 8)), 2: list(range(8, 15)), 3: list(range(15, 23)), 4: list(range(23, 28)), 5: list(range(28, 34)), 6: list(range(34, 41))}
OUT = {1: 'part1-foundations.html', 2: 'part2-signal-processing.html', 3: 'part3-physics-hardware.html',
       4: 'part4-image-formation.html', 5: 'part5-modes-capstones.html'}
TITLE = {1: 'Part I · Math & Signal Foundations', 2: 'Part II · Ultrasound Signal Processing & Estimation',
         3: 'Part III · Physics & Hardware', 4: 'Part IV · Image Formation & Quality', 5: 'Part V · Imaging Modes & Capstones'}
arg = sys.argv[1] if len(sys.argv) > 1 else 'b1'
if arg == 'all': arg = 'b1'
if arg in ('b1', 'b2', 'b3', 'b4', 'b5'):
    batch = int(arg[1]); mods, name = BATCH[batch], OUT[batch]
else:
    n = int(arg[1:]); batch = [b for b in BATCH if n in BATCH[b]][0]; mods = [n]; name = f'{arg}-only.html'
shell = (root / 'shell.html').read_text()
rail = (root / f'rail{batch}.html').read_text() if (root / f'rail{batch}.html').exists() else ''
home = (root / f'home{batch}.html').read_text() if (root / f'home{batch}.html').exists() else ''
secs, scripts = [], []
for n in mods:
    f = root / f'm{n}.html'
    if not f.exists():
        continue
    t = f.read_text()
    sc = re.findall(r'<script\b[^>]*>.*?</script>', t, flags=re.S)
    body = re.sub(r'<script\b[^>]*>.*?</script>', '', t, flags=re.S)
    secs.append(f'<!-- ===== {f.name} ===== -->\n' + body.strip())
    scripts.extend(sc)
import html as _html
def inline_code(t):
    def rep(m):
        p = root / m.group(1)
        src = p.read_text() if p.exists() else f'% MISSING FILE: {m.group(1)}'
        return f'<code data-file="{m.group(1)}">' + _html.escape(src.rstrip() + '\n', quote=False) + '</code>'
    return re.sub(r'<code data-file="([^"]+)">\s*</code>', rep, t)
secs = [inline_code(x) for x in secs]
out = (shell.replace('<!--RAIL-->', rail).replace('<!--HOME-->', home)
       .replace('<!--MODULES-->', '\n\n'.join(secs)).replace('<!--MODULE-SCRIPTS-->', '\n'.join(scripts)))
LOCAL = '--local' in sys.argv
_mj = (root / 'vendor' / 'mathjax-tex-svg-3.2.2.js')
if LOCAL and _mj.exists():
    out = out.replace('<!--MATHJAX-->', '<script>' + _mj.read_text().replace('</script', '<\\/script') + '</script>')
else:
    out = out.replace('<!--MATHJAX-->', '<script async src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>')
import json as _json
_u = _json.loads((root / 'urls.json').read_text()) if (root / 'urls.json').exists() else {}
for k, v in _u.items(): out = out.replace('{{' + k + '}}', v)
if batch in TITLE:
    out = out.replace('<title>Signal Processing Lab</title>', '<title>' + TITLE[batch] + '</title>', 1)
    out = out.replace("document.title = 'Signal Processing Lab';", "document.title = '" + TITLE[batch] + "';")
    out = out.replace("'sp-lab-view'", "'sp-lab%d-view'" % batch)
(root / 'build').mkdir(exist_ok=True)
p = root / 'build' / name
p.write_text(out)
print(p, f'{len(out)/1024:.0f} KB', f'{len(secs)} fragment(s)')
