#!/usr/bin/env python3
"""Headless test: Demo 25.6 default and tuned runs, live reduced ADMIRE timing in both modes, trade-off figure screenshot.
usage: pw7c.py [width] [scheme] [--nolive]"""
import sys, time, pathlib
from playwright.sync_api import sync_playwright
root = pathlib.Path('./'); out = pathlib.Path(__file__).parent / 'shots'; out.mkdir(exist_ok=True)
width = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 1280
scheme = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] in ('light', 'dark') else 'light'
page_file = root / 'build' / '_wrapped_m25-only.html'
errs = []
ROWS = "() => [...document.querySelectorAll('#m25-g-tbody tr')].map(r => r.innerText.replace(/\\s+/g, ' ')).join('\\n')"
def wait_done(pg):
    pg.wait_for_function("() => { const e = document.getElementById('m25-g-prog'); return e && e.textContent.startsWith('done'); }", timeout=300000)
def wait_sweeps(pg):
    pg.wait_for_function("() => document.getElementById('m25-a-prog').textContent.startsWith('swept') && document.getElementById('m25-s-prog').textContent.startsWith('swept')", timeout=300000)
def live(pg, tag):
    pg.wait_for_function("() => !document.getElementById('m25-g-adm-run').disabled", timeout=120000)
    print(tag, 'button label:', pg.inner_text('#m25-g-adm-run'))
    t1 = time.time(); pg.click('#m25-g-adm-run')
    pg.wait_for_function("() => document.getElementById('m25-g-adm-o').textContent.includes('computed live')", timeout=900000)
    print(f'{tag} live run wall {time.time() - t1:.1f} s;', pg.inner_text('#m25-g-adm-o'))
    print(tag, 'ADMIRE row:', pg.evaluate("() => [...document.querySelectorAll('#m25-g-tbody tr')].map(r => r.innerText.replace(/\\s+/g, ' ')).filter(t => t.startsWith('ADMIRE')).join(' || ')"))
with sync_playwright() as p:
    b = p.chromium.launch(); ctx = b.new_context(viewport={'width': width, 'height': 1000}, color_scheme=scheme, device_scale_factor=1)
    pg = ctx.new_page(); pg.route('**/fonts.g*/**', lambda r: r.abort())
    pg.on('console', lambda m: errs.append(f'console.{m.type}: {m.text}') if m.type in ('error', 'warning') and 'Failed to load resource' not in m.text else None)
    pg.on('pageerror', lambda e: errs.append(f'pageerror: {e}'))
    t0 = time.time(); pg.goto(page_file.as_uri() + '#m25'); wait_done(pg)
    print(f'default run done after {time.time() - t0:.1f} s:', pg.inner_text('#m25-g-prog')[:160])
    print('set readout:', pg.inner_text('#m25-g-set-o'))
    print(pg.evaluate(ROWS))
    el = pg.query_selector('#m25-t-cv').evaluate_handle("e => e.closest('.demo')").as_element(); el.scroll_into_view_if_needed(); pg.wait_for_timeout(500)
    el.screenshot(path=str(out / f'tune-{width}-{scheme}.png'))
    pg.select_option('#m25-t-y', 'full'); pg.wait_for_timeout(200); el.screenshot(path=str(out / f'tune-full-{width}-{scheme}.png')); pg.select_option('#m25-t-y', 'zoom')
    for v in ('adaptive', 'weights', 'other', 'all'): pg.select_option('#m25-t-set', v); pg.wait_for_timeout(100)
    pg.query_selector('#m25-protocol').screenshot(path=str(out / f'protocol-{width}-{scheme}.png'))
    wait_sweeps(pg)
    if '--nolive' not in sys.argv: live(pg, 'default')
    pg.select_option('#m25-g-set', 'tuned'); print('set readout after select:', pg.inner_text('#m25-g-set-o'))
    t0 = time.time(); pg.click('#m25-g-run'); pg.wait_for_timeout(300); wait_done(pg)
    print(f'tuned run done after {time.time() - t0:.1f} s:', pg.inner_text('#m25-g-prog')[:160])
    print('set readout:', pg.inner_text('#m25-g-set-o'))
    print(pg.evaluate(ROWS))
    g = pg.query_selector('#m25-g-das').evaluate_handle("e => e.closest('.demo')").as_element(); g.scroll_into_view_if_needed(); pg.wait_for_timeout(500); g.screenshot(path=str(out / f'g-tuned-{width}-{scheme}.png'))
    wait_sweeps(pg)
    if '--nolive' not in sys.argv: live(pg, 'tuned')
    print('errors:', errs if errs else 'none')
