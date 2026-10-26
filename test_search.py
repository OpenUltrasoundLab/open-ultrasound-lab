#!/usr/bin/env python3
"""Playwright test for the site search (needs a build in public/). Usage: python3 test_search.py [--shots]
Types four terms, asserts >=1 result each and that clicking the first result lands on a page containing the term;
also checks keyboard use (/, arrows, Enter, Esc), Greek tolerance, and takes screenshots when --shots is given."""
import sys, pathlib, threading, http.server, socketserver, functools
from playwright.sync_api import sync_playwright
root = pathlib.Path(__file__).resolve().parent
PUB = root / 'public'
TERMS = ['grating lobe', 'Kasai', 'CRLB', 'matching layer']

def serve():
    class Q(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *a): pass
    srv = socketserver.ThreadingTCPServer(('127.0.0.1', 0), functools.partial(Q, directory=str(PUB)))
    srv.daemon_threads = True
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv

def main():
    shots = '--shots' in sys.argv
    srv = serve(); base = f'http://127.0.0.1:{srv.server_address[1]}/'
    bad = []
    def check(ok, msg):
        print(('[ok  ] ' if ok else '[FAIL] ') + msg)
        if not ok: bad.append(msg)
    with sync_playwright() as p:
        b = p.chromium.launch()
        for width, scheme in [(1280, 'light'), (400, 'dark')]:
            ctx = b.new_context(viewport={'width': width, 'height': 900}, color_scheme=scheme)
            pg = ctx.new_page()
            errs = []
            pg.on('pageerror', lambda e: errs.append(str(e)))
            pg.route('**/*', lambda r: r.abort() if r.request.url.startswith('http') and '127.0.0.1' not in r.request.url else r.continue_())
            for term in TERMS:
                pg.goto(base + 'index.html')
                box = pg.locator('#os-q-hero')
                box.click(); box.fill(term)
                pg.wait_for_selector('#os-list a[role=option]', timeout=5000)
                n = pg.locator('#os-list a[role=option]').count()
                check(n >= 1, f'{width}px {scheme}: "{term}" -> {n} results')
                if shots and term == 'grating lobe':
                    pg.screenshot(path=str(root / 'shots' / f'search-{width}-{scheme}.png'))
                    pg.locator('#os-q-top').scroll_into_view_if_needed()
                pg.keyboard.press('Enter')   # first result is pre-selected
                pg.wait_for_load_state('load'); pg.wait_for_timeout(300)
                body = pg.evaluate('document.body.innerText').lower()
                hit = term.lower() in body or all(w in body for w in term.lower().split())
                check(hit and '/modules/' in pg.url, f'{width}px {scheme}: Enter on first "{term}" result -> {pg.url.replace(base, "")} contains term')
            # click navigation, keyboard, greek
            pg.goto(base + 'index.html'); box = pg.locator('#os-q-hero'); box.fill('Kasai')
            pg.wait_for_selector('#os-list a[role=option]')
            pg.locator('#os-list a[role=option]').first.click(); pg.wait_for_load_state('load')
            check('/modules/' in pg.url and 'kasai' in pg.evaluate('document.body.innerText').lower(), f'{width}px {scheme}: click on first Kasai result navigates')
            pg.goto(base + 'interview/index.html')
            pg.keyboard.press('/')
            check(pg.evaluate('document.activeElement.id') == 'os-q-top', 'interview page: "/" focuses the header search')
            pg.keyboard.type('theta'); pg.wait_for_selector('#os-list a[role=option]')
            n1 = pg.locator('#os-list a[role=option]').count()
            pg.keyboard.press('ArrowDown'); pg.keyboard.press('ArrowDown')
            check(pg.locator('#os-list a[aria-selected=true]').count() == 1 and pg.evaluate("document.activeElement.getAttribute('aria-activedescendant')") == 'os-opt-2', 'ArrowDown moves the selection')
            pg.keyboard.press('Escape')
            check(pg.locator('#os-list').is_hidden(), 'Esc closes the list')
            pg.fill('#os-q-top', 'θ'); pg.wait_for_selector('#os-list a[role=option]')
            check(pg.locator('#os-list a[role=option]').count() >= 1, f'Greek "θ" matches ("theta": {n1} results)')
            if shots and width == 400:
                pg.fill('#os-q-top', 'grating lobe'); pg.wait_for_timeout(300); pg.screenshot(path=str(root / 'shots' / 'search-400-dark-interview.png'))
            # module page: search box present (rail on desktop, meta row on phone), no horizontal overflow with list open
            pg.goto(base + 'modules/23-' + next(d.name[3:] for d in (PUB / 'modules').iterdir() if d.name.startswith('23-')) + '/index.html')
            vis = pg.locator('.os-q:visible').count()
            check(vis >= 1, f'{width}px module page has a visible search box')
            pg.locator('.os-q:visible').first.fill('grating lobe'); pg.wait_for_selector('#os-list a[role=option]')
            check(pg.evaluate('document.documentElement.scrollWidth - innerWidth') <= 0, f'{width}px no overflow with list open')
            if shots: pg.screenshot(path=str(root / 'shots' / f'search-module-{width}-{scheme}.png'))
            check(not errs, f'{width}px {scheme}: no JS errors {errs[:2]}')
            ctx.close()
        b.close()
    srv.shutdown()
    print(f'{len(bad)} failing'); sys.exit(1 if bad else 0)

if __name__ == '__main__':
    main()
