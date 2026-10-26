#!/usr/bin/env python3
"""Playwright render check over every HTML file in public/ (1280 px light, 400 px dark).
Fails on: JS/console errors, horizontal overflow, blank canvases, mjx-merror, raw TeX in visible text, MathJax missing.
Needs a build that works offline:  python3 build_site.py --local   (CDN builds are checked in CI, which has network access;
pass --cdn to skip the MathJax-presence check when the CDN is reachable).
Usage: python3 check_site.py [--fast] [--jobs N] [--only substring]   (--fast is accepted; the two viewport combos are always the minimum)
Screenshots: shots/site/<page>-<width>-<scheme>.png"""
import sys, pathlib, threading, http.server, socketserver, functools, argparse
from concurrent.futures import ProcessPoolExecutor
root = pathlib.Path(__file__).resolve().parent
PUB = root / 'public'
COMBOS = [(1280, 'light'), (400, 'dark')]

JS = """async () => {
  try { if (window.MathJax && MathJax.startup && MathJax.startup.promise) await MathJax.startup.promise; } catch (e) {}
  await new Promise(r => setTimeout(r, 900));
  const over = document.documentElement.scrollWidth - window.innerWidth;
  const skip = 'mjx-assistive-mml,mjx-container > svg,.eq,.eqn-math,.tbl-wrap,.readout,.pathline,.code pre,.mx-scroll,pre.cite,.rail,.iv-toc,.keyeq,.kmod';
  const wide = [];
  document.querySelectorAll('main *, body > * ').forEach(el => { const r = el.getBoundingClientRect(); if (r.width && r.right > window.innerWidth + 1 && !el.closest(skip)) wide.push(el.tagName.toLowerCase() + (el.id ? '#' + el.id : '') + (typeof el.className === 'string' && el.className ? '.' + el.className.split(' ')[0] : '')); });
  const blank = [];
  document.querySelectorAll('canvas').forEach(cv => { if (!cv.offsetParent) return; try { const g = cv.getContext('2d'); const d = g.getImageData(0, 0, cv.width, cv.height).data; let n = 0; for (let i = 3; i < d.length; i += 4 * 97) if (d[i]) { n++; if (n > 3) break; } if (n <= 3) blank.push(cv.id || '(no id)'); } catch (e) {} });
  const clone = document.body.cloneNode(true); clone.querySelectorAll('pre,code,script,style,mjx-container,details:not([open]) > :not(summary)').forEach(e => e.remove());
  const txt = clone.textContent;
  const raw = ['\\\\(', '\\\\[', '\\\\frac', '\\\\sum', '$$'].reduce((n, s) => n + txt.split(s).length - 1, 0);
  return { over, wide: wide.slice(0, 6), blank, merr: document.querySelectorAll('mjx-merror').length, raw,
           canvases: document.querySelectorAll('canvas').length, mj: !!(window.MathJax && MathJax.typesetPromise), needsMj: /\\\\\\(|\\\\\\[/.test(document.documentElement.innerHTML) };
}"""

def serve():
    h = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(PUB))
    h.log_message = lambda *a, **k: None
    class Q(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *a): pass
    srv = socketserver.ThreadingTCPServer(('127.0.0.1', 0), functools.partial(Q, directory=str(PUB)))
    srv.daemon_threads = True
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv

def run(args):
    files, cdn = args
    from playwright.sync_api import sync_playwright
    srv = serve(); base = f'http://127.0.0.1:{srv.server_address[1]}/'
    shots = root / 'shots' / 'site'; shots.mkdir(parents=True, exist_ok=True)
    out = []
    with sync_playwright() as p:
        b = p.chromium.launch()
        for width, scheme in COMBOS:
            ctx = b.new_context(viewport={'width': width, 'height': 900}, color_scheme=scheme)
            pg = ctx.new_page()
            def blocked(r):
                u = r.request.url
                return r.abort() if (u.startswith('http') and '127.0.0.1' not in u and not (cdn and 'jsdelivr' in u)) else r.continue_()
            pg.route('**/*', blocked)
            errs = []
            pg.on('console', lambda m: errs.append('console.error: ' + m.text) if m.type == 'error' and 'Failed to load resource' not in m.text and 'ERR_FAILED' not in m.text else None)
            pg.on('pageerror', lambda e: errs.append(f'pageerror: {e}'))
            for f in files:
                errs.clear()
                pg.goto(base + f)
                pg.wait_for_timeout(1200)
                info = pg.evaluate(JS)
                name = f.replace('/index.html', '').replace('.html', '').replace('/', '_') or 'index'
                try:
                    pg.screenshot(path=str(shots / f'{name}-{width}-{scheme}.png'), full_page=True, timeout=60000)
                except Exception as e:
                    print('       (screenshot skipped:', name, str(e)[:60], ')')
                probs = list(errs[:5])
                if info['over'] > 0: probs.append(f"overflow {info['over']}px")
                if info['wide']: probs.append(f"wide {info['wide']}")
                if info['blank']: probs.append(f"blank canvases {info['blank']}")
                if info['merr']: probs.append(f"mjx-merror x{info['merr']}")
                if info['raw']: probs.append(f"raw TeX x{info['raw']}")
                if info['needsMj'] and not info['mj'] and not cdn: probs.append('MathJax not loaded (build with --local)')
                out.append((f, width, scheme, info['canvases'], probs))
            ctx.close()
        b.close()
    srv.shutdown()
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--fast', action='store_true'); ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--only', default=''); ap.add_argument('--cdn', action='store_true')
    a = ap.parse_args()
    files = sorted(str(f.relative_to(PUB)) for f in PUB.rglob('*.html') if a.only in str(f))
    if not files: sys.exit('no html files in public/ (run build_site.py first)')
    chunks = [files[i::a.jobs] for i in range(a.jobs) if files[i::a.jobs]]
    bad = 0
    with ProcessPoolExecutor(len(chunks)) as ex:
        for res in ex.map(run, [(c, a.cdn) for c in chunks]):
            for f, w, sc, nc, probs in res:
                bad += bool(probs)
                print(f"[{'FAIL' if probs else 'ok  '}] {f:55s} {w}px {sc:5s} canvases={nc}" + ('  ' + '; '.join(probs) if probs else ''))
    print(f'{len(files)} files x {len(COMBOS)} viewports: {bad} failing')
    sys.exit(1 if bad else 0)
