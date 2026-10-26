#!/usr/bin/env python3
"""admire_grid.json + admire_imgs.json -> PNG sprite (8 noise columns x 3 aberration rows of 81 x 161 tiles, 8-bit grey)
-> <script type="application/json" id="m25-adm-data"> embedded in m25.html (replaced if present)."""
import json, sys, base64, io, re, pathlib
from PIL import Image
W = pathlib.Path(__file__).parent
meta = json.load(open(W / 'admire_grid.json')); imgs = json.load(open(W / 'admire_imgs.json'))
ms_per_px = float(sys.argv[1]) if len(sys.argv) > 1 else None
nx, nz, nc, nr = meta['nx'], meta['nz'], len(meta['noise']), len(meta['ab'])
sprite = Image.new('L', (nc * nx, nr * nz), 0)
for s, im in zip(meta['settings'], imgs):
    ci, ri = meta['noise'].index(s['nd']), meta['ab'].index(s['ab'])
    tile = Image.frombytes('L', (nx, nz), bytes(im))
    sprite.paste(tile, (ci * nx, ri * nz))
buf = io.BytesIO(); sprite.save(buf, format='PNG', optimize=True)
png = buf.getvalue()
# verify lossless round trip
back = Image.open(io.BytesIO(png)); assert back.mode == 'L'
for s, im in zip(meta['settings'], imgs):
    ci, ri = meta['noise'].index(s['nd']), meta['ab'].index(s['ab'])
    assert list(back.crop((ci * nx, ri * nz, ci * nx + nx, ri * nz + nz)).getdata()) == im
meta['png'] = 'data:image/png;base64,' + base64.b64encode(png).decode()
if ms_per_px: meta['msPerPixel'] = ms_per_px
blob = json.dumps(meta, separators=(',', ':'))
assert '</' not in blob
tag = '<script type="application/json" id="m25-adm-data">' + blob + '</script>\n'
p = pathlib.Path('./m25.html'); s = p.read_text()
if 'id="m25-adm-data"' in s:
    s = re.sub(r'<script type="application/json" id="m25-adm-data">.*?</script>\n', lambda m: tag, s, count=1, flags=re.S)
else:
    i = s.index('<script>\n(function () {\n/* ===== Module 25 engine')
    s = s[:i] + tag + s[i:]
p.write_text(s)
print(f'PNG {len(png)/1e3:.1f} kB, blob {len(blob)/1e3:.1f} kB, complete={meta["complete"]}')
