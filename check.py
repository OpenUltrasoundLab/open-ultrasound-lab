#!/usr/bin/env python3
"""Render-check the lab in headless Chromium.
Usage: python3 check.py m3            -> builds build/m3-only.html, checks view #m3
       python3 check.py b1 | all      -> checks build/part1-foundations.html (Part I), every view
       python3 check.py b2 .. b5      -> checks Part II .. Part V, every view
       add --fast to run only 1280-light and 400-dark
Reports console errors, page errors, horizontal overflow, and canvases that stayed blank.
Writes screenshots to shots/<view>-<width>-<scheme>.png (full page).
"""
import sys, subprocess, pathlib, json
from playwright.sync_api import sync_playwright
root = pathlib.Path(__file__).parent
arg = sys.argv[1] if len(sys.argv) > 1 else 'all'
if arg in ('all', 'b1', 'b2', 'b3', 'b4', 'b5'):
    b = 'b1' if arg == 'all' else arg
    subprocess.run([sys.executable, str(root / 'assemble.py'), b, '--local'], check=True)
    page_file = root / 'build' / {'b1': 'part1-foundations.html', 'b2': 'part2-signal-processing.html', 'b3': 'part3-physics-hardware.html', 'b4': 'part4-image-formation.html', 'b5': 'part5-modes-capstones.html'}[b]
    rng = {'b1': range(1, 8), 'b2': range(8, 15), 'b3': range(15, 23), 'b4': range(23, 28), 'b5': range(28, 34)}[b]
    views = ['home'] + [f'm{n}' for n in rng if (root / f'm{n}.html').exists()]
else:
    subprocess.run([sys.executable, str(root / 'assemble.py'), arg, '--local'], check=True)
    page_file = root / 'build' / f'{arg}-only.html'
    views = [arg]
shots = root / 'shots'; shots.mkdir(exist_ok=True)
html = '<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover"></head><body>' + page_file.read_text() + '</body></html>'
wrapped = root / 'build' / ('_wrapped_' + page_file.name)
wrapped.write_text(html)
problems = 0
with sync_playwright() as p:
    b = p.chromium.launch()
    combos = [(1280, 'light'), (400, 'dark')] if '--fast' in sys.argv else [(w, sc) for w in (1280, 400) for sc in ('light', 'dark')]
    for width, scheme in combos:
        if True:
            ctx = b.new_context(viewport={'width': width, 'height': 900}, color_scheme=scheme, device_scale_factor=1)
            pg = ctx.new_page()
            pg.route('**/fonts.g*/**', lambda r: r.abort())
            errs = []
            pg.on('console', lambda m: errs.append(f'console.{m.type}: {m.text}') if m.type in ('error',) and 'Failed to load resource' not in m.text else None)
            pg.on('pageerror', lambda e: errs.append(f'pageerror: {e}'))
            for v in views:
                errs.clear()
                pg.goto(wrapped.as_uri() + '#' + v)
                pg.wait_for_timeout(2500)
                info = pg.evaluate("""(v) => {
                  const sec = document.getElementById(v);
                  const over = document.documentElement.scrollWidth - window.innerWidth;
                  const wide = [];
                  if (sec) sec.querySelectorAll('*').forEach(el => { const r = el.getBoundingClientRect(); if (r.width && r.right > window.innerWidth + 1 && !el.closest('.eq,.eqn-math,.tbl-wrap,.readout,.pathline,.code pre,.mx-scroll')) wide.push(el.tagName.toLowerCase() + (el.id ? '#' + el.id : '') + (el.className && typeof el.className==='string' ? '.' + el.className.split(' ')[0] : '')); });
                  const blank = [];
                  if (sec) sec.querySelectorAll('canvas').forEach(cv => { try { const g = cv.getContext('2d'); const d = g.getImageData(0, 0, cv.width, cv.height).data; let n = 0; for (let i = 3; i < d.length; i += 4 * 97) if (d[i]) { n++; if (n > 3) break; } if (n <= 3) blank.push(cv.id || '(no id)'); } catch (e) {} });
                  const merr = sec ? sec.querySelectorAll('mjx-merror').length : 0; const raw = sec ? ['\\\\(', '\\\\[', '\\\\frac', '\\\\sum'].reduce((n, s) => n + sec.innerText.split(s).length - 1, 0) : 0;
                  return { visible: !!sec && !sec.hidden, over, wide: wide.slice(0, 8), blank, merr, raw, canvases: sec ? sec.querySelectorAll('canvas').length : 0 };
                }""", v)
                pg.screenshot(path=str(shots / f'{("home-" + page_file.stem) if v == "home" else v}-{width}-{scheme}.png'), full_page=True)
                bad = errs or info['over'] > 0 or info['wide'] or info['blank'] or not info['visible'] or info['merr'] or info['raw']
                problems += bool(bad)
                print(f"[{'FAIL' if bad else 'ok  '}] {v:5s} {width}px {scheme:5s} canvases={info['canvases']} overflow={info['over']}px"
                      + (f" wide={info['wide']}" if info['wide'] else '') + (f" blank={info['blank']}" if info['blank'] else '')
                      + (f" visible={info['visible']}" if not info['visible'] else '') + (f" texerr={info['merr']}" if info['merr'] else '') + (f" rawtex={info['raw']}" if info['raw'] else ''))
                for e in errs[:10]:
                    print('       ', e[:300])
            ctx.close()
    b.close()
print('screenshots in', shots)
sys.exit(1 if problems else 0)
