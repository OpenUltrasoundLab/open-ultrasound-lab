#!/usr/bin/env python3
"""Build the formula sheet: build/formula-sheet.html
Usage: python3 build_sheet.py            -> MathJax from the jsDelivr CDN (default)
       python3 build_sheet.py --local    -> MathJax embedded from vendor/ (same as assemble.py --local)

Pieces: sheet_body.html (tables, part structure), sheet_vis.html (convolution demo, rapid-fire cards),
the shell's <head> material and first <style> block, and the key equations of every module, read from
the .eqn.key boxes of m1.html ... m33.html at build time (marker <!--KEYEQ:a-b--> in sheet_body.html).
"""
import re, sys, json, pathlib
root = pathlib.Path(__file__).parent
LOCAL = '--local' in sys.argv

b = (root / 'sheet_body.html').read_text(encoding="utf-8")
v = (root / 'sheet_vis.html').read_text(encoding="utf-8")
conv = v[v.index('<!--CONV-->') + 11:v.index('<!--/CONV-->')]
rapid = v[v.index('<!--RAPID-->') + 12:v.index('<!--/RAPID-->')]
tail = v[v.index('<!--/RAPID-->') + 13:]

PART = lambda n: 1 if n <= 7 else 2 if n <= 14 else 3 if n <= 22 else 4 if n <= 27 else 5


def balanced_div(s, start):
    """Return the end index of the <div ...> element that starts at s[start]."""
    depth, i = 0, start
    for m in re.finditer(r'<div\b|</div>', s[start:]):
        depth += 1 if m.group(0) == '<div' else -1
        if depth == 0:
            return start + m.end()
    raise ValueError('unbalanced div at %d' % start)


def module_keyeqs(n):
    s = (root / f'm{n}.html').read_text(encoding="utf-8")
    title = re.search(r'data-title="([^"]*)"', s).group(1)
    out = []
    for m in re.finditer(r'<div class="eqn key"', s):
        blk = s[m.start():balanced_div(s, m.start())]
        tex = re.search(r'<div class="eqn-math">(.*?)</div>', blk, re.S).group(1).strip()
        no = re.search(r'class="eqn-no">\(([^)]*)\)', blk)
        note = re.search(r'<span class="note">(.*?)</span>\s*(?=</div>\s*$|<dl|<div)', blk, re.S)
        where = re.search(r'<dl class="where">.*?</dl>', blk, re.S)
        out.append((no.group(1) if no else '', tex, note.group(1).strip() if note else '', where.group(0) if where else ''))
    return title, out


def keyeq_html(a, z):
    h = []
    for n in range(a, z + 1):
        if not (root / f'm{n}.html').exists():
            continue
        title, eqs = module_keyeqs(n)
        h.append(f'  <h3 class="kmod" id="k-m{n}"><a href="{{{{PART{PART(n)}_URL}}}}" target="_blank" rel="noopener"><span class="n">{n:02d}</span>{title}</a></h3>')
        h.append('  <div class="keyeq">')
        for no, tex, note, where in eqs:
            # ids inside copied markup would collide with nothing here, but strip any to be safe
            note = re.sub(r'\sid="[^"]*"', '', note)
            h.append('    <div class="eqn key"><div class="eqn-body"><div class="eqn-math">' + tex + '</div>'
                     + (f'<span class="eqn-no">({no})</span>' if no else '') + '</div>'
                     + (f'<span class="note">{note}</span>' if note else '')
                     + (f'<details class="sym"><summary>symbols</summary>{where}</details>' if where else '')
                     + '</div>')
        h.append('  </div>')
    return '\n'.join(h)


counts = {}
def fill(m):
    a, z = int(m.group(1)), int(m.group(2))
    counts[(a, z)] = sum(len(module_keyeqs(n)[1]) for n in range(a, z + 1) if (root / f'm{n}.html').exists())
    return keyeq_html(a, z)
b = re.sub(r'<!--KEYEQ:(\d+)-(\d+)-->', fill, b)

# convolution demo after the theorems; rapid-fire cards replace the rapid-fire heading
b = b.replace('<!--CONV-->', conv)
i = b.index('<h2 id="s-rapid">'); j = b.index('<p class="muted">', i); b = b[:i] + rapid + '\n  ' + b[j:]

n = [0]
def num(m): n[0] += 1; return '<span class="sec">%d</span>' % n[0]
b = re.sub(r'<span class="sec">\d+</span>', num, b)

shell = (root / 'shell.html').read_text(encoding="utf-8")
head_bits = shell[shell.index('<link'):shell.index('</style>') + len('</style>')]
mj_cfg = re.search(r'<script>\s*window\.MathJax = .*?</script>', shell, re.S).group(0)
mj_cfg = re.sub(r'\nwindow\.Lab &&.*?\n', '\n', mj_cfg)
mj = root / 'vendor' / 'mathjax-tex-svg-3.2.2.js'
if LOCAL and mj.exists():
    mj_tag = '<script>' + mj.read_text(encoding="utf-8").replace('</script', '<\\/script') + '</script>'
else:
    mj_tag = '<script async src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>'

out = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
       '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n'
       '<title>Ultrasound Formula Sheet</title>\n' + head_bits + '\n</head>\n<body>\n'
       + b + '\n' + tail + '\n' + mj_cfg + '\n' + mj_tag + '\n</body>\n</html>\n')
urls = json.loads((root / 'urls.json').read_text(encoding="utf-8")) if (root / 'urls.json').exists() else {}
for k, val in urls.items():
    out = out.replace('{{' + k + '}}', val)
(root / 'build').mkdir(exist_ok=True)
p = root / 'build' / 'formula-sheet.html'
p.write_text(out, encoding="utf-8")
print(p, f'{len(out)/1024:.0f} KB', 'sections', n[0], 'key equations', sum(counts.values()), counts,
      'MathJax', 'embedded' if (LOCAL and mj.exists()) else 'CDN')
