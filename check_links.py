#!/usr/bin/env python3
"""Link checker for public/: every internal href/src must resolve to a file, and #fragments to an id/name in the target page.
Usage: python3 check_links.py [public_dir]   (exit 1 on any broken link)"""
import sys, pathlib, re
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path(__file__).parent / 'public').resolve()

class P(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True); self.links = []; self.ids = set()
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if a.get('id'): self.ids.add(a['id'])
        if tag == 'a' and a.get('name'): self.ids.add(a['name'])
        for k in ('href', 'src'):
            if a.get(k) is not None and tag != 'base': self.links.append((tag, k, a[k]))

cache = {}
def parse(p):
    if p not in cache:
        q = P(); q.feed(p.read_text(errors='replace')); cache[p] = q
    return cache[p]

bad, n = [], 0
for f in sorted(root.rglob('*.html')):
    for tag, k, v in parse(f).links:
        v = v.strip()
        if not v or re.match(r'^(https?:|mailto:|tel:|javascript:|data:|//)', v): continue
        n += 1
        u = urlsplit(v); path = unquote(u.path)
        if not path:
            tgt = f
        else:
            tgt = (root / path.lstrip('/')) if path.startswith('/') else (f.parent / path)
            tgt = pathlib.Path(posix := str(tgt)).resolve()
            if tgt.is_dir(): tgt = tgt / 'index.html'
            if not tgt.is_file() or root not in tgt.parents:
                bad.append(f'{f.relative_to(root)}: <{tag} {k}="{v}"> -> missing file'); continue
        if u.fragment and u.fragment != 'top' and tgt.suffix == '.html':
            if unquote(u.fragment) not in parse(tgt).ids:
                bad.append(f'{f.relative_to(root)}: <{tag} {k}="{v}"> -> no id "{u.fragment}" in {tgt.relative_to(root)}')
for b in bad[:200]: print('BROKEN', b)
print(f'checked {n} internal links in {len(cache) and len(list(root.rglob("*.html")))} html files: {len(bad)} broken')
sys.exit(1 if bad else 0)
