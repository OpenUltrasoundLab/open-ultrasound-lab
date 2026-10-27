#!/usr/bin/env python3
"""Build the public static site for Open Ultrasound Lab into public/.

Usage: python3 build_site.py            -> public/ with MathJax from the jsDelivr CDN (for GitHub Pages)
       python3 build_site.py --local    -> same, but MathJax served from public/vendor/ (offline testing)

Reads at build time: shell.html, m1.html ... m33.html, rail1-5.html (short nav titles), home1-5.html
(part ledes, card subtitles), matlab/*.m, AUDIT_*.md (private) or audit_summary.json, errata.json, CITATION.cff, site.config.json and
the output of build_sheet.py. Never edits any of those. Re-run after any module edit.

Layout of public/:
  index.html                      landing (parts map, learning paths, citation)
  modules/NN-slug/index.html      one page per module
  interview/index.html            every module's interview questions
  formula-sheet/index.html        build_sheet.py output
  about.html disclaimer.html errata.html cite.html contribute.html
  matlab/*.m  CITATION.cff  .nojekyll  (vendor/ only with --local)
"""
import re, sys, json, shutil, html, posixpath, pathlib, subprocess, unicodedata

root = pathlib.Path(__file__).resolve().parent
PUB = root / 'public'
LOCAL = '--local' in sys.argv
MJ_FILE = 'mathjax-tex-svg-3.2.2.js'
MJ_CDN = 'https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js'

# ---------------------------------------------------------------- config
CFG = json.loads((root / 'site.config.json').read_text())
def cfg_expand(s):
    for _ in range(3):
        for k, v in CFG.items():
            if isinstance(v, str):
                s = s.replace('{{' + k + '}}', v)
    return s
for k in list(CFG):
    if isinstance(CFG[k], str):
        CFG[k] = cfg_expand(CFG[k])
VERSION, DATE = CFG['VERSION'], CFG['DATE']
HAS_DOI = bool(CFG.get('DOI')) and '0000000' not in CFG['DOI']   # placeholder DOI = cite by URL until Zenodo mints one
ARCH = f"Zenodo. https://doi.org/{CFG['DOI']}" if HAS_DOI else CFG['COURSE_URL']
REPO, AUTHOR = CFG['REPO_URL'], CFG['AUTHOR_URL']

# ---------------------------------------------------------------- course structure
PARTS = [  # (roman, name, first, last)
    ('I', 'Math & signal foundations', 1, 7),
    ('II', 'Ultrasound signal processing & estimation', 8, 14),
    ('III', 'Physics & hardware', 15, 22),
    ('IV', 'Image formation & quality', 23, 27),
    ('V', 'Imaging modes & capstones', 28, 33),
]
def part_of(n):
    return next(i + 1 for i, p in enumerate(PARTS) if p[2] <= n <= p[3])

PATHS = [
    ('ulm', 'Ultrasound localization microscopy',
     'From phasors and SVD clutter filtering to microbubble detection, localization, tracking and resolution metrics.',
     [1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 23, 25, 26, 28, 31]),
    ('sos', 'Sound speed and aberration',
     'Delay estimation, inverse problems and optimization, then coherence, k-space and the two capstones.',
     [1, 2, 3, 4, 6, 8, 11, 12, 13, 14, 23, 24, 32, 33]),
    ('hw', 'Transducer and hardware engineer',
     'The signal basics, then wave physics, piezoelectrics, stacks, CMUT/PMUT, arrays, front ends, characterization and acoustic output.',
     [1, 2, 3, 4] + list(range(15, 23))),
    ('int', 'Interview preparation',
     'Every module’s interview and qualifying-exam questions with model answers, on one page, grouped by module.',
     None),
]


def esc(s):
    return html.escape(s, quote=False)


def attr(s):
    return html.escape(s, quote=True)


def slugify(t):
    t = unicodedata.normalize('NFKD', t).encode('ascii', 'ignore').decode()
    t = t.lower().replace('&', ' and ')
    words = re.findall(r'[a-z0-9]+', t)
    out = ''
    for w in words:
        if len(out) + len(w) + 1 > 60:
            break
        out = (out + '-' + w) if out else w
    return out


def plain(s):
    """HTML fragment -> plain one-line text (TeX delimiters dropped)."""
    s = re.sub(r'<[^>]+>', '', s)
    s = html.unescape(s)
    s = s.replace('\\(', '').replace('\\)', '')
    return re.sub(r'\s+', ' ', s).strip()


def clip(s, n=158):
    if len(s) <= n:
        return s
    s = s[:n].rsplit(' ', 1)[0].rstrip(',;:')
    return s + '…'


def balanced(s, start, tag):
    depth = 0
    for m in re.finditer(r'<%s\b|</%s>' % (tag, tag), s[start:]):
        depth += 1 if not m.group(0).startswith('</') else -1
        if depth == 0:
            return start + m.end()
    raise ValueError('unbalanced <%s> at %d' % (tag, start))


# ---------------------------------------------------------------- read modules
short = {}
for k in range(1, 6):
    f = root / f'rail{k}.html'
    if f.exists():
        for m in re.finditer(r'data-view="m(\d+)"><span class="num">\d+</span><span class="t">(.*?)</span>', f.read_text()):
            short[int(m.group(1))] = html.unescape(m.group(2))
subs, part_lede = {}, {}
for k in range(1, 6):
    f = root / f'home{k}.html'
    if not f.exists():
        continue
    t = f.read_text()
    for m in re.finditer(r'<p class="sub">(.*?)</p>.*?href="#m(\d+)"', t, re.S):
        subs[int(m.group(2))] = html.unescape(m.group(1))
    m = re.search(r'<p class="lede">(.*?)</p>', t, re.S)
    if m:
        part_lede[k] = plain(m.group(1))

MODS = {}
for n in range(1, 34):
    f = root / f'm{n}.html'
    if not f.exists():
        continue
    s = f.read_text()
    title = html.unescape(re.search(r'data-title="([^"]*)"', s).group(1))
    lede_m = re.search(r'<p class="lede">(.*?)</p>', s, re.S)
    MODS[n] = dict(n=n, src=s, title=title, short=short.get(n, title), sub=subs.get(n, ''),
                   lede=plain(lede_m.group(1)) if lede_m else '', slug=f'{n:02d}-{slugify(title)}')
    MODS[n]['url'] = f"modules/{MODS[n]['slug']}/index.html"
if not MODS:
    sys.exit('no module files found')


# ---------------------------------------------------------------- link rewriting
def rel(target, from_dir):
    """target and from_dir are paths relative to public/; returns a relative URL."""
    r = posixpath.relpath(target, from_dir or '.')
    return r


def rewrite_links(t, from_dir, cur):
    """Rewrite <a href> in a module fragment. cur = module number whose content this is (for #home etc.);
    from_dir = directory of the output page relative to public/. Page of module `cur` is 'self' only
    when from_dir equals that module's directory."""
    self_dir = posixpath.dirname(MODS[cur]['url'])
    stats = {'n': 0}

    def modhref(n, frag=''):
        if n not in MODS:
            return None
        if posixpath.dirname(MODS[n]['url']) == from_dir:
            return frag or '#top'
        return rel(MODS[n]['url'], from_dir) + frag

    def fix_tag(m):
        tag, inner = m.group(1), m.group(2)
        hm = re.search(r'\bhref="([^"]*)"', tag)
        if not hm:
            return m.group(0)
        h = hm.group(1)
        new = None
        if h == '#home':
            new = rel('index.html', from_dir) + f'#part-{part_of(cur)}'
        elif re.fullmatch(r'#m\d+', h):
            new = modhref(int(h[2:]))
        elif re.fullmatch(r'#m\d+-[\w-]+', h):
            n = int(re.match(r'#m(\d+)', h).group(1))
            new = modhref(n, h)
        elif h.startswith('{{PART'):
            pm = re.match(r'\{\{PART(\d)_URL\}\}(#.*)?$', h)
            txt = plain(inner)
            mm = re.search(r'\bModules? (\d{1,2})\b', txt) or re.match(r'^(?:←\s*)?Part [IVX]+ · (\d{1,2})\b', txt)
            frag = pm.group(2) or '' if pm else ''
            if mm and int(mm.group(1)) in MODS:
                n = int(mm.group(1))
                new = modhref(n, frag if frag.startswith(f'#m{n}-') else '')
            elif pm:
                new = rel('index.html', from_dir) + f'#part-{pm.group(1)}'
        if new is None:
            return m.group(0)
        stats['n'] += 1
        tag = tag.replace(hm.group(0), f'href="{new}"')
        tag = re.sub(r'\s+target="_blank"', '', tag)
        tag = re.sub(r'\s+rel="noopener"', '', tag)
        return tag + inner

    t = re.sub(r'(<a\b[^>]*>)(.*?)(?=</a>)', fix_tag, t, flags=re.S)
    return t, stats['n']


def inline_code(t):
    def rep(m):
        p = root / m.group(1)
        src = p.read_text() if p.exists() else f'% MISSING FILE: {m.group(1)}'
        return f'<code data-file="{m.group(1)}">' + esc(src.rstrip() + '\n') + '</code>'
    return re.sub(r'<code data-file="([^"]+)">\s*</code>', rep, t)


# ---------------------------------------------------------------- shared chrome
shell = (root / 'shell.html').read_text()
HEAD_CSS = shell[shell.index('<link'):shell.index('</style>') + len('</style>')]
MJ_CFG = re.search(r'<script>\s*window\.MathJax = .*?</script>', shell, re.S).group(0)


def mj_tag(depth):
    if LOCAL:
        return f'<script async src="{"../" * depth}vendor/{MJ_FILE}"></script>'
    return f'<script async src="{MJ_CDN}"></script>'


SITE_CSS = r'''<style>
/* ---- Open Ultrasound Lab: site chrome (added by build_site.py; tokens from the shell) ---- */
.skip { position: absolute; left: -999px; } .skip:focus { left: 12px; top: 12px; z-index: 50; background: var(--panel); padding: 6px 10px; }
.site-rail .brand { text-transform: none; font-stretch: 82%; font-size: 1.12rem; }
.site-rail .part h4 { line-height: 1.3; }
.site-rail details.part > summary { cursor: pointer; list-style: none; font: 600 0.72rem/1.3 var(--f-display); letter-spacing: 0.12em; text-transform: uppercase; color: var(--muted); padding: 2px 0; }
.site-rail details.part > summary::-webkit-details-marker { display: none; }
.site-rail details.part > summary::before { content: "\25B8\00a0"; }
.site-rail details.part[open] > summary::before { content: "\25BE\00a0"; }
.site-rail details.part nav { margin-top: 6px; }
.site-rail .quick a { display: block; padding: 3px 8px; font: 500 0.88rem/1.4 var(--f-display); text-decoration: none; border-radius: var(--radius); }
.site-rail .quick a:hover { background: var(--grid); }
.site-rail .allmods { display: none; }
.page-meta { display: flex; justify-content: flex-end; gap: 16px; max-width: 980px; font: 500 0.8rem/1.3 var(--f-display); margin: -10px 0 8px; }
.page-meta a { color: var(--muted); }
.page-meta a:hover { color: var(--accent); }
.dl-bar { font: 400 0.86rem/1.6 var(--f-display); color: var(--muted); background: var(--panel); border: 1px solid var(--rule); border-radius: var(--radius); padding: 8px 12px; margin: 10px 0; }
.dl-bar a { font-family: var(--f-mono); font-size: 0.82rem; white-space: nowrap; }
.site-foot { max-width: 980px; margin-top: 64px; padding-top: 16px; border-top: 1px solid var(--rule); font: 400 0.84rem/1.6 var(--f-display); color: var(--muted); display: flex; flex-wrap: wrap; gap: 4px 18px; }
.site-foot a { color: var(--muted); } .site-foot a:hover { color: var(--accent); }
.site-foot .ver { color: var(--ink); font-weight: 600; }
@media (max-width: 860px) {
  .site-rail details.part, .site-rail .quick, .site-rail .part h4 { display: none; }
  .site-rail .part.cur { display: contents; }
  .site-rail .allmods { display: block; white-space: nowrap; font: 600 0.82rem var(--f-display); text-decoration: none; padding: 6px 8px; }
  .page-meta { margin-top: 0; }
}
/* ---- standalone pages (landing, about, ...) ---- */
.top { display: flex; flex-wrap: wrap; align-items: baseline; gap: 6px 22px; padding: 14px clamp(16px, 4vw, 56px); border-bottom: 1px solid var(--rule); }
.top .brand { font: 700 1.08rem/1.15 var(--f-display); font-stretch: 82%; color: var(--ink); text-decoration: none; margin-right: auto; }
.top nav { display: flex; flex-wrap: wrap; gap: 4px 18px; font: 500 0.9rem/1.4 var(--f-display); }
.top nav a { color: var(--ink); text-decoration: none; } .top nav a:hover, .top nav a[aria-current="page"] { color: var(--accent); }
.page { max-width: 1060px; padding: 28px clamp(16px, 4vw, 56px) 64px; }
.page .prose { max-width: 72ch; }
.page pre.cite { font: 400 0.8rem/1.55 var(--f-mono); background: var(--code-bg); border: 1px solid var(--rule); border-radius: var(--radius); padding: 12px 14px; overflow-x: auto; white-space: pre; }
.page .site-foot { max-width: none; }
.errata-entry { background: var(--panel); border: 1px solid var(--rule); border-radius: 8px; padding: 12px 16px; margin: 14px 0; }
.errata-entry h3 { margin: 0 0 4px; }
.flag-list { font-size: 0.92rem; } .flag-list li { margin: 6px 0; max-width: 80ch; }
.flag-list .where { font: 600 0.78rem var(--f-mono); color: var(--muted); }
/* landing */
.hero { padding: 36px 0 6px; max-width: 980px; }
.hero h1 { font-size: clamp(2.6rem, 8vw, 5.2rem); font-stretch: 68%; line-height: 0.95; margin: 0 0 14px; }
.hero .lede { max-width: 58ch; }
.hero .free { font: 600 1.05rem/1.4 var(--f-display); color: var(--accent); margin: 14px 0 4px; }
.hero .by { color: var(--muted); font-size: 0.95rem; margin: 2px 0 0; }
.hero .cta { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 18px; }
.hero .cta a { text-decoration: none; }
.aperture-wrap { margin: 34px 0 6px; }
.aperture-cap { font: 400 0.84rem/1.45 var(--f-display); color: var(--muted); margin: 0 0 8px; max-width: 72ch; }
.aperture { display: grid; grid-template-columns: repeat(33, minmax(0, 1fr)); gap: 3px; }
.aperture a.el { display: flex; align-items: flex-end; justify-content: center; height: 78px; border-radius: 3px 3px 1px 1px; text-decoration: none;
  background: color-mix(in srgb, var(--pc) 24%, var(--panel)); color: var(--ink); font: 600 0.7rem/1 var(--f-mono); padding-bottom: 6px; transition: background .25s, color .25s, transform .25s; }
.aperture a.el:hover, .aperture a.el:focus-visible { background: var(--pc); color: var(--panel); }
:root { --part1: #1F5FD1; --part2: #0F7A63; --part3: #A86E05; --part4: #C2366F; --part5: #6B48C8; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --part1: #7FA3FF; --part2: #4CCFAE; --part3: #E3B04B; --part4: #F07DAE; --part5: #B39CFF; } }
:root[data-theme="dark"] { --part1: #7FA3FF; --part2: #4CCFAE; --part3: #E3B04B; --part4: #F07DAE; --part5: #B39CFF; }
.p1 { --pc: var(--part1); } .p2 { --pc: var(--part2); } .p3 { --pc: var(--part3); } .p4 { --pc: var(--part4); } .p5 { --pc: var(--part5); }
.tier { border-left: 4px solid var(--pc, var(--accent)); padding: 2px 12px 6px 14px; background: color-mix(in srgb, var(--pc, var(--accent)) 6%, transparent); border-radius: 0 10px 10px 0; }
.tier h3 { color: var(--pc, inherit); }
.tier li .num, .site-rail .part .num { color: var(--pc, var(--muted)) !important; }
.site-rail .part h4, .site-rail .part summary { border-left: 3px solid var(--pc, transparent); padding-left: 7px; }
.hero .eyebrow { font: 600 0.85rem/1.2 var(--f-display); letter-spacing: 0.12em; text-transform: uppercase; color: var(--accent); margin: 0 0 10px; }
.hero .stats { list-style: none; display: flex; flex-wrap: wrap; gap: 10px; padding: 0; margin: 18px 0 16px; }
.hero .stats li { display: flex; flex-direction: column; min-width: 112px; padding: 10px 14px; border: 1px solid var(--rule); border-radius: 10px; background: var(--panel); }
.hero .stats b { font: 700 1.55rem/1.1 var(--f-display); font-variant-numeric: tabular-nums; }
.hero .stats span { color: var(--muted); font-size: 0.85rem; }
.note { position: relative; display: grid; grid-template-columns: 56px minmax(0, 1fr); gap: 4px 12px; max-width: 820px; margin: 32px 0 10px; padding: 22px 24px 20px 14px; border-radius: 14px; background: var(--panel); border: 1px solid var(--rule); border-left: 6px solid transparent; border-image: linear-gradient(180deg, var(--part1), var(--part2), var(--part3), var(--part4), var(--part5)) 1; }
.note-mark { font: 700 4.4rem/0.9 Georgia, serif; color: var(--accent); opacity: 0.35; text-align: center; }
.note-body h2 { margin: 0 0 6px; font-size: 1.3rem; }
.note-lead { font: 600 1.15rem/1.45 var(--f-display); margin: 0 0 8px; }
.note-body p { margin: 0 0 10px; max-width: 66ch; }
.note-points { list-style: none; display: flex; flex-wrap: wrap; gap: 8px; padding: 0; margin: 4px 0 12px; }
.note-points li { padding: 6px 12px; border-radius: 999px; background: color-mix(in srgb, var(--accent) 10%, transparent); font-size: 0.92rem; }
.note-sig { display: flex; flex-wrap: wrap; align-items: center; gap: 6px 12px; margin-top: 6px; padding-top: 10px; border-top: 1px dashed var(--rule); }
.note-sig .name { font: italic 600 1.1rem/1.2 Georgia, serif; }
.note-sig .role { color: var(--muted); font-size: 0.9rem; margin-right: auto; }
.btn.small { padding: 4px 12px; font-size: 0.85rem; text-decoration: none; }
@media (max-width: 560px) { .note { grid-template-columns: minmax(0, 1fr); padding: 18px 16px; } .note-mark { display: none; } .hero .stats { display: grid; grid-template-columns: 1fr 1fr; } .hero .stats li { min-width: 0; } }
.why { max-width: 760px; margin: 28px 0 6px; padding: 16px 20px; border: 1px solid var(--rule); border-radius: 12px; background: var(--panel); }
.why h2 { margin: 0 0 6px; font-size: 1.2rem; } .why p { margin: 0 0 8px; } .why .sig { color: var(--muted); font-style: italic; margin: 0; }
.aperture.has-path a.el { background: var(--grid); color: var(--muted); }
.aperture.has-path a.el.on { background: var(--pc); color: var(--panel); transform: translateY(-6px); }
.aperture .plab { grid-row: 2; font: 600 0.74rem/1.25 var(--f-display); color: var(--muted); border-top: 2px solid var(--pc); padding-top: 5px; min-width: 0; }

.pathbar { display: flex; flex-wrap: wrap; gap: 8px; margin: 16px 0 4px; }
.pathbar .btn[aria-pressed="true"] { background: var(--accent); color: var(--panel); }
.paths { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 14px; margin-top: 14px; }
.path { background: var(--panel); border: 1px solid var(--rule); border-radius: 10px; padding: 14px 16px; min-width: 0; }
.path h3 { margin: 0 0 4px; }
.path p { margin: 0 0 8px; font-size: 0.93rem; }
.path ol { margin: 0; padding-left: 0; list-style: none; font-size: 0.9rem; line-height: 1.5; columns: 2 12em; column-gap: 18px; }
.path ol li { break-inside: avoid; }
.path ol .num { font: 600 0.76rem var(--f-mono); color: var(--muted); margin-right: 6px; }
.paths.js .path[hidden] { display: none; }
.tier li .num { font: 600 0.78rem/1.6 var(--f-mono); color: var(--muted); }
.tier li small { display: block; color: var(--muted); font-size: 0.82rem; }
.diff { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 14px 28px; margin-top: 10px; }
.diff > div { min-width: 0; } .diff h3 { margin-top: 6px; } .diff p { margin: 0 0 8px; font-size: 0.95rem; }
.diff .prov { margin: 0 4px 4px 0; }
.iv-mod { border-top: 1px solid var(--rule); margin-top: 30px; }
.iv-mod h2 { font-size: var(--step-2); margin-top: 18px; }
.iv-mod h2 a { color: inherit; text-decoration-color: var(--accent); }
.iv-mod h2 .n { font: 600 0.9rem/1 var(--f-mono); color: var(--muted); margin-right: 10px; vertical-align: 0.2em; }
.iv-toc { columns: 3 15em; column-gap: 24px; font-size: 0.9rem; padding-left: 0; list-style: none; }
.iv-toc .num { font: 600 0.76rem var(--f-mono); color: var(--muted); margin-right: 6px; }
@media (max-width: 700px) {
  .aperture { grid-template-columns: repeat(11, minmax(0, 1fr)); row-gap: 6px; }
  .aperture a.el { height: 46px; }
  .aperture .plab { display: none; }
}
@media (prefers-reduced-motion: reduce) { .aperture a.el { transition: none; } .aperture.has-path a.el.on { transform: none; } }
</style>'''


# ---------------------------------------------------------------- site search (index built from the module sources)
SEARCH_CSS = r"""
/* ---- site search (search.js) ---- */
.vh { position: absolute !important; width: 1px; height: 1px; overflow: hidden; clip: rect(0 0 0 0); white-space: nowrap; }
.os-box { position: relative; display: block; min-width: 0; }
.os-box input { width: 100%; box-sizing: border-box; font: 500 0.9rem/1.3 var(--f-display); color: var(--ink); background: var(--panel); border: 1px solid var(--rule);
  border-radius: var(--radius); padding: 7px 10px; -webkit-appearance: none; appearance: none; }
.os-box input::placeholder { color: var(--muted); opacity: 1; }
.os-box input:focus-visible, .os-box input:focus { outline: 2px solid var(--accent); outline-offset: 1px; border-color: var(--accent); }
.top .os-box { width: 15rem; }
.os-rail { margin: 8px 0 10px; }
.os-meta { display: none; margin-right: auto; flex: 1 1 12rem; max-width: 22rem; }
.os-hero { max-width: 34rem; margin: 18px 0 0; }
.os-hero input { font-size: 1.05rem; padding: 12px 14px; border-width: 2px; }
.os-kbd { position: absolute; right: 8px; top: 50%; transform: translateY(-50%); font: 600 0.72rem/1 var(--f-mono); color: var(--muted); border: 1px solid var(--rule); border-radius: 4px; padding: 2px 6px; pointer-events: none; }
.os-box input:focus + .os-kbd, .os-box.has-q .os-kbd { display: none; }
#os-list { position: fixed; z-index: 100; overflow-y: auto; overscroll-behavior: contain; background: var(--panel); color: var(--ink); border: 1px solid var(--rule);
  border-radius: 8px; box-shadow: 0 10px 30px rgba(0,0,0,.22); padding: 4px 0; font: 400 0.88rem/1.4 var(--f-display); }
#os-list[hidden] { display: none; }
#os-list .os-grp { padding: 8px 12px 3px; font: 700 0.74rem/1.3 var(--f-display); letter-spacing: .06em; text-transform: uppercase; color: var(--muted); }
#os-list .os-grp b { font-family: var(--f-mono); margin-right: 6px; color: var(--accent); }
#os-list a { display: block; padding: 6px 12px; color: var(--ink); text-decoration: none; border-left: 3px solid transparent; }
#os-list a[aria-selected="true"], #os-list a:hover { background: var(--grid); border-left-color: var(--accent); }
#os-list a:focus-visible { outline: 2px solid var(--accent); outline-offset: -2px; }
#os-list .os-t { display: inline-block; font: 700 0.66rem/1 var(--f-mono); text-transform: uppercase; letter-spacing: .04em; color: var(--panel); background: var(--muted);
  border-radius: 3px; padding: 2px 5px; margin-right: 7px; vertical-align: 1px; }
#os-list .os-t.t-mod { background: var(--accent); } #os-list .os-t.t-sec { background: var(--sig4); } #os-list .os-t.t-demo { background: var(--sig2); }
#os-list .os-t.t-eq { background: var(--sig3); } #os-list .os-t.t-iv, #os-list .os-t.t-chk { background: var(--ink); }
#os-list mark { background: color-mix(in srgb, var(--sig3) 38%, transparent); color: inherit; border-radius: 2px; padding: 0 1px; }
#os-list .os-s { overflow-wrap: anywhere; }
#os-list .os-none { padding: 10px 12px; color: var(--muted); }
:target { scroll-margin-top: 76px; }
.demo:target, .rule:target, .eqn:target, details.quiz:target, tr:target { outline: 2px solid var(--accent); outline-offset: 2px; }
@media (max-width: 860px) { .site-rail .os-rail { display: none; } .os-meta { display: block; } .top .os-box { width: 100%; order: 3; } }
"""
SITE_CSS = SITE_CSS.replace('</style>', SEARCH_CSS + '</style>')


def search_box(kind, depth):
    ph = {'hero': 'Search the course: equations, demos, interview questions…'}.get(kind, 'Search (press /)')
    return (f'<div class="os-box os-{kind}" role="search"><label class="vh" for="os-q-{kind}">Search the course</label>'
            f'<input type="search" id="os-q-{kind}" class="os-q" placeholder="{attr(ph)}" autocomplete="off" spellcheck="false" '
            f'role="combobox" aria-expanded="false" aria-controls="os-list" aria-autocomplete="list" aria-haspopup="listbox">'
            + ('' if kind == 'hero' else '<kbd class="os-kbd" aria-hidden="true">/</kbd>') + '</div>')


def search_script(depth):
    up = '../' * depth
    return f'<script defer src="{up}search.js" data-root="{up}"></script>\n'


GREEK = {'alpha': 'α', 'beta': 'β', 'gamma': 'γ', 'delta': 'δ', 'epsilon': 'ε', 'varepsilon': 'ε', 'zeta': 'ζ', 'eta': 'η', 'theta': 'θ', 'vartheta': 'θ',
         'kappa': 'κ', 'lambda': 'λ', 'mu': 'μ', 'nu': 'ν', 'xi': 'ξ', 'pi': 'π', 'rho': 'ρ', 'sigma': 'σ', 'tau': 'τ', 'phi': 'φ', 'varphi': 'φ',
         'chi': 'χ', 'psi': 'ψ', 'omega': 'ω', 'Delta': 'Δ', 'Gamma': 'Γ', 'Sigma': 'Σ', 'Omega': 'Ω', 'Phi': 'Φ', 'Psi': 'Ψ', 'Lambda': 'Λ', 'Theta': 'Θ'}
SYM = {'le': '≤', 'leq': '≤', 'ge': '≥', 'geq': '≥', 'mp': '∓', 'pm': '±', 'circ': '°', 'cdot': '·', 'times': '×', 'approx': '≈', 'sim': '~', 'infty': '∞',
       'rightarrow': '→', 'Rightarrow': '⇒', 'to': '→', 'neq': '≠', 'ne': '≠', 'propto': '∝', 'sqrt': '√', 'sum': 'Σ', 'int': '∫', 'partial': '∂', 'nabla': '∇', 'll': '≪', 'gg': '≫'}
_TEX_DROP = r'mathbf|mathrm|mathcal|mathbb|boldsymbol|text|textbf|operatorname|left|right|hat|bar|tilde|widehat|vec|big|Big|displaystyle|lVert|rVert|quad|qquad'


def ptxt(s, limit=0):
    """HTML fragment -> plain text for the search index (TeX reduced to readable symbols)."""
    t = plain(s)
    t = t.replace('\\[', '').replace('\\]', '')
    if '\\' in t:
        t = re.sub(r'\\[lr]?[Vv]ert\b', '|', t)
        t = re.sub(r'\\(?:%s)(?![A-Za-z])' % _TEX_DROP, '', t)
        t = re.sub(r'\\([A-Za-z]+)', lambda m: GREEK.get(m.group(1), SYM.get(m.group(1), m.group(1))), t)
        t = re.sub(r'\\[,;:! ]', ' ', t)
    t = t.replace('{', '').replace('}', '').replace('\\', '')
    t = re.sub(r'\s+', ' ', t).strip()
    return t[:limit] if limit else t


SEARCH_INDEX = {'mods': {}, 'e': []}


# ---------- glossary (glossary.json): first use of a term per section gets a pop-up card ----------
GLOSS = json.loads((root / 'glossary.json').read_text()) if (root / 'glossary.json').exists() else []
GL_BY = {g['t']: g for g in GLOSS}
GL_RX = re.compile(r'(?<![\w\-/.#&])(' + '|'.join(re.escape(t) for t in sorted(GL_BY, key=len, reverse=True)) + r')(?=s?(?![\w\-]))') if GLOSS else None
GL_SKIP = {'script', 'style', 'pre', 'code', 'a', 'h1', 'h2', 'h3', 'h4', 'summary', 'canvas', 'svg', 'button', 'select', 'option', 'textarea', 'label', 'abbr', 'title'}
GL_FIRST = {}   # term -> (module, section id) of first tagged use in its teaching module


def gl_id(t):
    return 'g-' + re.sub(r'[^A-Za-z0-9]+', '-', t).strip('-')


def tag_terms(body, n):
    """Wrap the first use of each glossary term in every h2 section (text only; not math, code, links, headings)."""
    if not GL_RX:
        return body, set()
    out, used, seen, stack = [], set(), set(), []
    sec = 'top'
    in_math = False
    for m in re.finditer(r'(<[^>]+>)|([^<]+)', body):
        tag, txt = m.group(1), m.group(2)
        if tag:
            tm = re.match(r'<(/?)([a-zA-Z0-9]+)([^>]*?)(/?)>', tag)
            if tm:
                close, name, attrs, selfc = tm.group(1), tm.group(2).lower(), tm.group(3), tm.group(4)
                if name == 'h2' and not close:
                    seen = set(); hid = re.search(r'id="([^"]+)"', attrs); sec = hid.group(1) if hid else sec
                skip = name in GL_SKIP or re.search(r'class="[^"]*\b(eqn-math|where|tag|sec|num|page-meta|dl-bar)\b', attrs or '')
                if not close and not selfc and name not in ('br', 'img', 'input', 'hr', 'meta', 'link', 'wbr', 'source', 'col'):
                    stack.append((name, bool(skip) or (stack[-1][1] if stack else False)))
                elif close:
                    for k in range(len(stack) - 1, -1, -1):
                        if stack[k][0] == name:
                            del stack[k:]; break
            out.append(tag); continue
        if stack and stack[-1][1]:
            out.append(txt); continue
        # split off TeX: \( \) \[ \] $$ ; state carries across nodes
        parts = re.split(r'(\\\(|\\\)|\\\[|\\\]|\$\$)', txt)
        res = []
        for part in parts:
            if part in ('\\(', '\\['): in_math = True; res.append(part); continue
            if part in ('\\)', '\\]'): in_math = False; res.append(part); continue
            if part == '$$': in_math = not in_math; res.append(part); continue
            if in_math or not part:
                res.append(part); continue
            def rep(mm):
                t = mm.group(1); g = GL_BY[t]
                if t in seen or ('only' in g and n not in g['only']):
                    return t
                seen.add(t); used.add(t)
                if g['m'] == n and t not in GL_FIRST:
                    GL_FIRST[t] = (n, sec)
                return f'<abbr class="gl" tabindex="0" data-g="{t}">{t}</abbr>'
            res.append(GL_RX.sub(rep, part))
        out.append(''.join(res))
    return ''.join(out), used


def gl_link(t, depth):
    g = GL_BY[t]; up = '../' * depth
    n, sec = GL_FIRST.get(t, (g['m'], ''))
    return up + MODS[n]['url'].replace('index.html', '') + (('#' + sec) if sec else ''), n


GL_CSS = """<style>
abbr.gl { text-decoration: underline dotted; text-decoration-color: var(--muted); text-underline-offset: 3px; cursor: help; }
abbr.gl:hover, abbr.gl:focus { text-decoration-color: var(--accent); outline: none; }
.gl-card { position: absolute; z-index: 60; box-sizing: border-box; max-width: min(360px, calc(100vw - 32px)); background: var(--panel, #fff); color: var(--ink);
  border: 1px solid var(--rule); border-radius: 10px; box-shadow: 0 8px 28px rgba(0,0,0,.18); padding: 12px 14px; font: 400 0.9rem/1.45 var(--f-body); }
.gl-card b { font-family: var(--f-display); } .gl-card .f { color: var(--muted); }
.gl-card p { margin: 6px 0 8px; } .gl-card .l { display: flex; gap: 14px; font-size: 0.82rem; flex-wrap: wrap; }
.gl-card .x { position: absolute; top: 4px; right: 8px; border: 0; background: none; color: var(--muted); font-size: 1.1rem; cursor: pointer; }
.gl-list { list-style: none; padding: 0; margin: 0; } .gl-list li { padding: 10px 0; border-bottom: 1px solid var(--rule); scroll-margin-top: 80px; }
.gl-list li:target { background: var(--accent-wash); } .gl-list .t { font: 600 1rem var(--f-display); } .gl-list .f { color: var(--muted); }
.gl-list p { margin: 4px 0; } .gl-letters { display: flex; flex-wrap: wrap; gap: 6px 10px; margin: 8px 0 16px; }
#gl-filter { width: 100%; max-width: 420px; padding: 8px 10px; border: 1px solid var(--rule); border-radius: 8px; background: var(--paper); color: var(--ink); font: inherit; }
</style>"""

GL_JS = """<script>
(function () {
  var D = window.GL_DATA || {}, card = null, cur = null;
  function close() { if (card) { card.remove(); card = null; } if (cur) { cur.setAttribute('aria-expanded', 'false'); cur = null; } }
  function open(el) {
    var g = D[el.dataset.g]; if (!g) return; if (cur === el) return; close(); cur = el;
    card = document.createElement('div'); card.className = 'gl-card'; card.setAttribute('role', 'dialog');
    card.innerHTML = '<button class="x" aria-label="Close">\\u00d7</button><b></b> <span class="f"></span><p></p><div class="l"><a class="m"></a><a class="a">All terms</a></div>';
    card.querySelector('b').textContent = el.dataset.g; card.querySelector('.f').textContent = g[0];
    card.querySelector('p').textContent = g[1];
    var m = card.querySelector('.m'); m.href = g[2]; m.textContent = 'Taught in Module ' + g[3];
    card.querySelector('.a').href = g[4];
    card.querySelector('.x').onclick = close;
    document.body.appendChild(card); el.setAttribute('aria-expanded', 'true');
    var r = el.getBoundingClientRect(), w = card.offsetWidth, x = Math.max(window.scrollX + 16, Math.min(r.left + window.scrollX, window.scrollX + document.documentElement.clientWidth - w - 16));
    var y = r.bottom + window.scrollY + 6; if (r.bottom + card.offsetHeight + 12 > window.innerHeight && r.top > card.offsetHeight + 12) y = r.top + window.scrollY - card.offsetHeight - 6;
    card.style.left = x + 'px'; card.style.top = y + 'px';
  }
  document.addEventListener('click', function (e) { var a = e.target.closest && e.target.closest('abbr.gl'); if (a) { e.preventDefault(); cur === a ? close() : open(a); } else if (card && !card.contains(e.target)) close(); });
  document.addEventListener('keydown', function (e) { var a = e.target.closest && e.target.closest('abbr.gl'); if (a && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); open(a); } if (e.key === 'Escape') close(); });
  var hov; document.addEventListener('mouseover', function (e) { if (!window.matchMedia('(hover: hover)').matches) return; var a = e.target.closest && e.target.closest('abbr.gl'); if (a) { clearTimeout(hov); hov = setTimeout(function () { open(a); }, 250); } });
  document.addEventListener('mouseout', function (e) { var a = e.target.closest && e.target.closest('abbr.gl'); if (a) clearTimeout(hov); });
})();
</script>"""


def gl_data(terms, depth):
    d = {}
    for t in sorted(terms):
        g = GL_BY[t]; href, mn = gl_link(t, depth)
        d[t] = [g['f'], g['d'], href, mn, ('../' * depth) + 'glossary.html#' + gl_id(t)]
    return '<script>window.GL_DATA=' + json.dumps(d, ensure_ascii=False).replace('</', '<\\/') + ';</script>'


def add_anchors(body, n):
    """Give sections, demos, rules, equations, cheat-sheet rows and quiz questions stable ids (built output only)."""
    used = set(re.findall(r'\bid="([^"]+)"', body))
    def uniq(i):
        base, k = i, 2
        while i in used:
            i = f'{base}-{k}'; k += 1
        used.add(i); return i
    def h2(m):
        if m.group(1):
            return m.group(0)
        sec = re.search(r'<span class="sec">([\w.]+)</span>', m.group(2))
        i = uniq(f'm{n}-s' + (sec.group(1).replace('.', '-') if sec else slugify(plain(m.group(2)))))
        return f'<h2 id="{i}">{m.group(2)}</h2>'
    body = re.sub(r'<h2( id="[^"]*")?>(.*?)</h2>', h2, body, flags=re.S)
    def demo(m):
        tag = re.search(r'<span class="tag">Demo ([\d.]+)</span>', body[m.end():m.end() + 400])
        return f'<div class="demo" id="{uniq(f"m{n}-demo-" + (tag.group(1).replace(".", "-") if tag else str(len(used))))}"{m.group(1) or ""}>'
    body = re.sub(r'<div class="demo"(?! id=)([^>]*)>', demo, body)
    k = [0]
    def rule(m):
        k[0] += 1
        return f'<div class="rule" id="{uniq(f"m{n}-rule-{k[0]}")}">'
    body = re.sub(r'<div class="rule">', rule, body)
    def eqn(m):
        return f'<div class="{m.group(1)}" id="{uniq(f"m{n}-eqx{len(used)}")}">'
    body = re.sub(r'<div class="(eqn(?: key)?)">', eqn, body)
    def quiz(m):
        q = re.match(r'<summary><span class="q">([^<]+)</span>', body[m.end():])
        return f'<details class="quiz" id="{uniq(f"m{n}-q-" + (q.group(1) if q else str(len(used))))}"{m.group(1)}>'
    body = re.sub(r'<details class="quiz"(?! id=)([^>]*)>', quiz, body)
    # cheat-sheet rows
    for cm in list(re.finditer(r'<h2[^>]*>(?:(?!</h2>).)*Cheat sheet</h2>', body, re.S)):
        a = cm.end(); b = body.find('<h2', a); b = len(body) if b < 0 else b
        seg = body[a:b]; j = [0]
        def tr(m):
            j[0] += 1
            return f'<tr id="{uniq(f"m{n}-cs-{j[0]}")}">'
        body = body[:a] + re.sub(r'<tr>(?=<td)', tr, seg) + body[b:]
        break
    return body


def index_module(n, M, body):
    ents = SEARCH_INDEX['e']
    SEARCH_INDEX['mods'][str(n)] = [M['title'], M['url']]
    ents.append(['mod', ptxt(M['title'] + ' — ' + M['lede']), n, ''])
    heads = [(m.start(), m.group(1), ptxt(m.group(2))) for m in re.finditer(r'<h2 id="([^"]+)">(.*?)</h2>', body, re.S)]
    for _, i, t in heads:
        ents.append(['sec', re.sub(r'^(\d+\.\d+)(?=\S)', r'\1 ', t), n, i])
    for m in re.finditer(r'<div class="demo" id="([^"]+)">.{0,600}?<h3>(.*?)</h3>', body, re.S):
        tag = re.search(r'<span class="tag">(.*?)</span>', m.group(0))
        ents.append(['demo', ptxt((tag.group(1) + ': ' if tag else '') + m.group(2)), n, m.group(1)])
    for m in re.finditer(r'<div class="eqn( key)?" id="([^"]+)">', body):
        seg = body[m.start():balanced(body, m.start(), 'div')]
        wh = re.search(r'<dl class="where">(.*?)</dl>', seg, re.S)
        if not (m.group(1) or wh):
            continue
        tex = re.search(r'<div class="eqn-math">(.*?)</div>', seg, re.S)
        no = re.search(r'<span class="eqn-no">(.*?)</span>', seg)
        note = re.search(r'<span class="note">(.*?)</span>\s*$', seg.rsplit('</div>', 1)[0] + '</span>', re.S) or re.search(r'<span class="note">(.*?)</span>', seg, re.S)
        wtxt = ''
        if wh:
            wtxt = '; '.join(f'{ptxt(a)}: {ptxt(b)}' for a, b in re.findall(r'<dt>(.*?)</dt><dd>(.*?)</dd>', wh.group(1), re.S))
        main = ptxt(note.group(1)) if note else wtxt
        extra = ptxt(tex.group(1), 260) if tex else ''
        if note and wtxt:
            extra += ' | ' + wtxt
        ents.append(['eq', ((ptxt(no.group(1)) + ' ') if no else '') + main[:420], n, m.group(2), extra[:600]])
    for m in re.finditer(r'<div class="rule" id="([^"]+)">\s*<div class="say">(.*?)</div>', body, re.S):
        ents.append(['rule', ptxt(m.group(2)), n, m.group(1)])
    for m in re.finditer(r'<tr id="([^"]+-cs-\d+)">(.*?)</tr>', body, re.S):
        cells = [ptxt(c) for c in re.findall(r'<td[^>]*>(.*?)</td>', m.group(2), re.S)]
        if cells:
            ents.append(['cheat', ' · '.join(cells)[:300], n, m.group(1)])
    for m in re.finditer(r'<details class="quiz" id="([^"]+)"[^>]*><summary><span class="q">([^<]+)</span><span>(.*?)</span></summary>', body, re.S):
        sec = [t for pos, i, t in heads if pos < m.start()]
        kind = 'iv' if sec and 'Interview' in sec[-1] else 'chk'
        q = re.sub(r'<span class="muted">.*?</span>', '', m.group(3), flags=re.S)
        ents.append([kind, ptxt(q), n, m.group(1)])


def write_search_index():
    out = json.dumps({'v': 1, 'mods': SEARCH_INDEX['mods'], 'e': SEARCH_INDEX['e']}, ensure_ascii=False, separators=(',', ':'))
    (PUB / 'search-index.json').write_text(out, encoding='utf-8')
    shutil.copy2(root / 'search.js', PUB / 'search.js')
    return len(out.encode('utf-8')), len(SEARCH_INDEX['e'])


def footer(depth):
    up = '../' * depth
    return (f'<footer class="site-foot" role="contentinfo">'
            f'<span class="ver">Open Ultrasound Lab v{esc(VERSION)} · last reviewed {esc(DATE)}</span>'
            f'<span>Text and figures <a href="https://creativecommons.org/licenses/by/4.0/" rel="license">CC BY 4.0</a> · code MIT</span>'
            f'<a href="{up}glossary.html">Glossary</a><a href="{up}disclaimer.html">Disclaimer</a><a href="{up}errata.html">Errata</a>'
            f'<a href="{up}cite.html">Cite</a><a href="{up}about.html">About</a>'
            f'<a href="{up}contribute.html">Contribute</a>'
            f'<span>By M. Nguyen (<a href="{attr(AUTHOR)}" rel="noopener">LinkedIn</a>), with Claude (Anthropic)</span>'
            f'</footer>')


def head(title, desc, depth, extra=''):
    return ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n'
            f'<title>{esc(title)}</title>\n<meta name="description" content="{attr(desc)}">\n'
            f'<meta property="og:title" content="{attr(title)}">\n<meta property="og:description" content="{attr(desc)}">\n'
            + HEAD_CSS + '\n' + SITE_CSS + '\n' + search_script(depth) + extra + '</head>\n')


def topbar(depth, cur=''):
    up = '../' * depth
    items = [('index.html#map', 'Modules', 'map'), ('index.html#paths', 'Learning paths', 'paths'),
             ('formula-sheet/index.html', 'Formula sheet', 'sheet'), ('interview/index.html', 'Interview prep', 'interview'),
             ('glossary.html', 'Glossary', 'glossary'), ('about.html', 'About', 'about'), ('errata.html', 'Errata', 'errata')]
    nav = ''.join(f'<a href="{up}{h}"' + (' aria-current="page"' if k == cur else '') + f'>{t}</a>' for h, t, k in items)
    return (f'<a class="skip" href="#main">Skip to content</a><header class="top"><a class="brand" href="{up}index.html">Open Ultrasound Lab</a>'
            f'<nav aria-label="Site">{nav}</nav>{search_box("top", depth)}</header>\n')


def static_page(path, title, desc, body, cur='', mathjax=False):
    depth = path.count('/')
    extra = ''
    tail = ''
    if mathjax:
        cfg = re.sub(r'\nwindow\.Lab &&.*?\n', '\n', MJ_CFG)
        tail = cfg + '\n' + mj_tag(depth) + '\n'
    out = (head(title, desc, depth, extra) + '<body>\n' + topbar(depth, cur)
           + f'<main id="main" class="page">\n{body}\n{footer(depth)}\n</main>\n' + tail + '</body>\n</html>\n')
    write(path, out)


written = []
def write(path, text):
    p = PUB / path
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    written.append(path)


# ---------------------------------------------------------------- module pages
def rail_html(cur):
    from_dir = posixpath.dirname(MODS[cur]['url'])
    up = rel('.', from_dir) if from_dir else '.'
    up = up + '/'
    h = [f'<aside class="rail site-rail" aria-label="Course modules">',
         f'<a class="brand" href="{up}index.html">Open Ultrasound Lab<small>33 modules in five parts</small></a>',
         f'<div class="quick"><a href="{up}index.html#paths">Learning paths</a><a href="{up}formula-sheet/index.html">Formula sheet</a>'
         f'<a href="{up}interview/index.html">Interview prep</a></div>', search_box('rail', 2)]
    for i, (ro, name, a, z) in enumerate(PARTS, 1):
        links = []
        for n in range(a, z + 1):
            if n not in MODS:
                continue
            href = '#top' if n == cur else rel(MODS[n]['url'], from_dir)
            ac = ' aria-current="page"' if n == cur else ''
            links.append(f'<a class="nav" href="{href}" data-view="m{n}"{ac}><span class="num">{n:02d}</span><span class="t">{esc(MODS[n]["short"])}</span></a>')
        nav = '<nav>' + ''.join(links) + '</nav>'
        label = f'Part {ro} · {esc(name)}'
        if part_of(cur) == i:
            h.append(f'<div class="part cur p{i}"><h4>{label}</h4>{nav}</div>')
        else:
            h.append(f'<details class="part p{i}"><summary>{label}</summary>{nav}</details>')
    h.append(f'<a class="allmods" href="{up}index.html#map">All parts</a>')
    h.append('</aside>')
    return '\n'.join(h)


SHELL_FIXES = [
    ("function show(id) {", "const SINGLE = (document.body && document.body.dataset.single) || '';\n  function show(id) {"),
    ("try { localStorage.setItem('sp-lab-view', id); } catch (e) {}", "if (!SINGLE) { try { localStorage.setItem('sp-lab-view', id); } catch (e) {} }"),
    ("document.title = 'Signal Processing Lab';", ""),
    ("window.addEventListener('hashchange', () => { show(location.hash.slice(1) || 'home'); window.scrollTo(0, 0); });",
     "if (!SINGLE) window.addEventListener('hashchange', () => { show(location.hash.slice(1) || 'home'); window.scrollTo(0, 0); });"),
    ("let start = location.hash.slice(1); if (!start) {", "let start = SINGLE || location.hash.slice(1); if (!start) {"),
]
# keep a section anchor in view after fonts / MathJax change the layout
ANCHOR_JS = '''<script>
(function () { if (!location.hash || location.hash === '#top') return;
  var go = function () { try { var el = document.getElementById(decodeURIComponent(location.hash.slice(1))); if (el) el.scrollIntoView(); } catch (e) {} };
  window.addEventListener('load', function () { go(); setTimeout(go, 600);
    var t = setInterval(function () { if (window.MathJax && MathJax.startup && MathJax.startup.promise) { clearInterval(t); MathJax.startup.promise.then(go); } }, 100);
    setTimeout(function () { clearInterval(t); }, 8000); });
})();
</script>'''

shell_body = shell[shell.index('<div class="app">'):]
for a, b in SHELL_FIXES:
    if a not in shell_body:
        sys.exit('shell.html changed; cannot apply fix: ' + a[:60])
    shell_body = shell_body.replace(a, b, 1)

for _n, _M in MODS.items():   # glossary pre-pass: where each term is first used in its teaching module
    tag_terms(add_anchors(re.sub(r'<script\b[^>]*>.*?</script>', '', _M['src'], flags=re.S), _n), _n)

link_count = 0
for n, M in MODS.items():
    from_dir = posixpath.dirname(M['url'])
    s = M['src']
    scripts = re.findall(r'<script\b[^>]*>.*?</script>', s, flags=re.S)
    body = re.sub(r'<script\b[^>]*>.*?</script>', '', s, flags=re.S).strip()
    # single-module page: section visible without a hash
    body = re.sub(r'(<section class="view" id="m%d"[^>]*?)\s+hidden(\s*>)' % n, r'\1\2', body, count=1)
    # stable anchor on the interview heading
    body = re.sub(r'<h2>(<span class="sec">[\d.]+</span>Interview and qualifying-exam questions)</h2>',
                  f'<h2 id="m{n}-interview">\\1</h2>', body, count=1)
    body = add_anchors(body, n)
    index_module(n, M, body)
    body = inline_code(body)
    body, gl_used = tag_terms(body, n)
    files = re.findall(r'<code data-file="(matlab/[^"]+)"', body)
    if files:
        seen = list(dict.fromkeys(files))
        links = ' · '.join(f'<a href="../../{f}" download>{esc(f.split("/")[-1])}</a>' for f in seen)
        bar = (f'<p class="dl-bar">Download the scripts: {links}. Base MATLAB R2019b+ or GNU Octave 8, no toolboxes; '
               f'<a href="{attr(REPO)}/tree/main/matlab">all scripts in the repository</a>.</p>')
        i = body.index('<details class="code"')
        body = body[:i] + bar + '\n' + body[i:]
    body, c1 = rewrite_links(body, from_dir, n)
    scripts_txt, c2 = rewrite_links('\n'.join(scripts), from_dir, n)
    link_count += c1 + c2
    title = f'{n:02d} · {M["title"]} — Open Ultrasound Lab'
    issue = f'{REPO}/issues/new?title=Module%20{n:02d}:%20'
    meta = (f'<div class="page-meta" id="top">{search_box("meta", 2)}<a href="{attr(issue)}">Report an error</a>'
            f'<a href="../../errata.html">Errata</a></div>')
    main = meta + '\n' + body + '\n' + footer(2)
    out = shell_body.replace('<!--RAIL-->', rail_html(n)).replace('<!--HOME-->', '').replace('<!--MODULES-->', main)
    out = out.replace('<!--MODULE-SCRIPTS-->', scripts_txt + '\n' + ANCHOR_JS + '\n' + gl_data(gl_used, 2) + GL_JS)
    out = out.replace('<!--MATHJAX-->', mj_tag(2))
    out = head(title, clip(M['lede']), 2, GL_CSS) + f'<body data-single="m{n}">\n' + out + '\n</body>\n</html>\n'
    leftover = re.findall(r'\{\{[A-Z][A-Z0-9_]+\}\}', out)
    if leftover:
        sys.exit(f'm{n}: unresolved placeholders {sorted(set(leftover))}')
    write(M['url'], out)


# ---------------------------------------------------------------- interview page
def interview_page():
    parts_html, toc, total = [], [], 0
    for n, M in MODS.items():
        s = M['src']
        i = s.find('Interview and qualifying-exam questions')
        if i < 0:
            continue
        j = s.find('<h2', i + 10)
        seg = s[i:j if j > 0 else len(s)]
        qs = []
        for m in re.finditer(r'<details class="quiz"', seg):
            qs.append(seg[m.start():balanced(seg, m.start(), 'details')])
        if not qs:
            continue
        total += len(qs)
        qhtml, _ = rewrite_links('\n'.join(qs), 'interview', n)
        href = '../' + M['url'] + f'#m{n}-interview'
        toc.append(f'<li><span class="num">{n:02d}</span><a href="#iv-m{n}">{esc(M["short"])}</a></li>')
        parts_html.append(f'<section class="iv-mod" id="iv-m{n}"><h2><span class="n">{n:02d}</span>'
                          f'<a href="../{M["url"]}">{esc(M["title"])}</a></h2>'
                          f'<p class="muted">{len(qs)} questions · <a href="{href}">open in the module</a></p>\n{qhtml}</section>')
    body = (f'<div class="hero"><h1>Interview preparation</h1><p class="lede">{total} interview and qualifying-exam questions '
            f'from all {len(MODS)} modules, each with a model answer and the follow-up an interviewer would ask next. '
            'Tags mark the level: screen, onsite or research. Click a question to reveal the answer.</p></div>'
            f'<h2>Modules</h2><ul class="iv-toc">{"".join(toc)}</ul>\n' + '\n'.join(parts_html))
    static_page('interview/index.html', 'Interview preparation — Open Ultrasound Lab',
                f'{total} ultrasound interview and qualifying-exam questions with model answers, grouped by module.',
                body, 'interview', mathjax=True)
    return total


n_questions = interview_page()


# ---------------------------------------------------------------- landing page
def count(pat):
    return sum(len(re.findall(pat, M['src'])) for M in MODS.values())


n_demos = count(r'<div class="demo"')
n_eq = count(r'<div class="eqn(?: key)?" id=')
n_lit = count(r'<div class="litfig"')
n_scripts = len(list((root / 'matlab').glob('m*.m')))




def landing():
    els, labs = [], []
    for i, (ro, name, a, z) in enumerate(PARTS, 1):
        for n in range(a, z + 1):
            if n in MODS:
                ttl = attr(f'{n:02d} {MODS[n]["title"]}')
                els.append(f'<a class="el p{i}" href="{MODS[n]["url"]}" data-n="{n}" title="{ttl}" aria-label="{ttl}">{n:02d}</a>')
        labs.append(f'<span class="plab p{i}" style="grid-column: {a} / {z + 1}">Part {ro} · {esc(name)}</span>')
    aperture = ('<div class="aperture" id="aperture">' + ''.join(els) + ''.join(labs) + '</div>')

    tiers = []
    for i, (ro, name, a, z) in enumerate(PARTS, 1):
        items = ''.join(
            f'<li><span class="num">{n:02d}</span><span><a href="{MODS[n]["url"]}">{esc(MODS[n]["title"])}</a>'
            + (f'<small>{esc(MODS[n]["sub"])}</small>' if MODS[n]['sub'] else '') + '</span></li>'
            for n in range(a, z + 1) if n in MODS)
        lede = part_lede.get(i, '')
        first = clip(lede.split('. ')[0] + '.', 220) if lede else ''
        tiers.append(f'<div class="tier p{i}" id="part-{i}"><h3>Part {ro} · {esc(name)}</h3>'
                     f'<p class="sub">Modules {a:02d}–{z:02d}</p>'
                     + (f'<p class="muted" style="font-size:.9rem;margin:0 0 10px">{esc(first)}</p>' if first else '')
                     + f'<ol>{items}</ol></div>')

    paths, btns = [], []
    for key, name, desc, mods in PATHS:
        if mods is None:
            body = (f'<p>{esc(desc)}</p><p><a class="btn solid" href="interview/index.html" style="text-decoration:none;display:inline-block">'
                    f'Open the {n_questions} questions</a></p>')
            data = ''
        else:
            lis = ''.join(f'<li><span class="num">{n:02d}</span><a href="{MODS[n]["url"]}">{esc(MODS[n]["short"])}</a></li>' for n in mods if n in MODS)
            body = f'<p>{esc(desc)} {len(mods)} modules, in order.</p><ol>{lis}</ol>'
            data = ' data-mods="' + ','.join(str(n) for n in mods) + '"'
        paths.append(f'<div class="path" id="path-{key}"{data}><h3>{esc(name)}</h3>{body}</div>')
        btns.append(f'<button class="btn" type="button" data-path="{key}" aria-pressed="false" aria-controls="path-{key}">{esc(name)}</button>')

    js = '''<script>
(function () {
  var ap = document.getElementById('aperture'), wrap = document.getElementById('paths-list'), bar = document.getElementById('pathbar');
  if (!ap || !wrap || !bar) return;
  var cards = wrap.querySelectorAll('.path'), btns = bar.querySelectorAll('button');
  bar.hidden = false; wrap.classList.add('js');
  function pick(key) {
    var mods = [];
    btns.forEach(function (b) { var on = b.dataset.path === key; b.setAttribute('aria-pressed', on ? 'true' : 'false'); });
    cards.forEach(function (c) { var on = !key || c.id === 'path-' + key; c.hidden = !on; if (key && on && c.dataset.mods) mods = c.dataset.mods.split(',').map(Number); });
    ap.classList.toggle('has-path', mods.length > 0);
    ap.querySelectorAll('a.el').forEach(function (a) { a.classList.toggle('on', mods.indexOf(+a.dataset.n) >= 0); });
  }
  btns.forEach(function (b) { b.addEventListener('click', function () { pick(b.getAttribute('aria-pressed') === 'true' ? '' : b.dataset.path); }); });
  pick('');
})();
</script>'''

    body = f'''<div class="hero">
  <p class="eyebrow">Free · open · interactive</p>
  <h1>Open Ultrasound Lab</h1>
  <p class="lede">The mathematics, physics and signal processing behind modern ultrasound imaging, from phasors and the FFT to transducers, beamforming, localization microscopy, sound-speed estimation and aberration correction. Every idea comes with a figure you can change.</p>
  <ul class="stats" aria-label="Course at a glance">
    <li><b>{len(MODS)}</b><span>modules in 5 parts</span></li>
    <li><b>{n_demos}</b><span>live demos</span></li>
    <li><b>{n_eq}</b><span>numbered equations</span></li>
    <li><b>{n_scripts}</b><span>MATLAB/Octave labs</span></li>
    <li><b>{n_questions}</b><span>interview questions</span></li>
  </ul>
  {search_box('hero', 0)}
  <div class="cta"><a class="btn solid" href="{MODS[min(MODS)]["url"]}">Start with module 01</a><a class="btn" href="#paths">Choose a learning path</a><a class="btn" href="formula-sheet/index.html">Formula sheet</a></div>
  <p class="by">Version {esc(VERSION)} · {esc(DATE)}</p>
</div>
<section class="note" aria-labelledby="note-h">
  <div class="note-mark" aria-hidden="true">&ldquo;</div>
  <div class="note-body">
    <h2 id="note-h">Why I built this</h2>
    <p class="note-lead">I am an ultrasound engineer, and I am still learning every day.</p>
    <p>For years I wished the fundamentals, from phasors and the FFT to transducers, beamforming and the newest imaging methods, were gathered in one place and explained with figures you can change. This is my attempt, built with the help of AI and shared freely, especially for students and young engineers who are just getting started.</p>
    <ul class="note-points">
      <li><b>Free.</b> Not a product; nothing to buy.</li>
      <li><b>Checked.</b> Numbers come from code you can run.</li>
      <li><b>Imperfect.</b> It certainly contains mistakes; please tell me.</li>
    </ul>
    <p>If you can explain a topic better, please contribute. I hope it saves you some of the time it would have saved me.</p>
    <div class="note-sig"><span class="name">M. Nguyen</span><span class="role">Ultrasound engineer</span>
      <a class="btn small" href="{attr(AUTHOR)}" rel="noopener">LinkedIn</a><a class="btn small" href="{attr(REPO)}">GitHub</a></div>
  </div>
</section>
<div class="aperture-wrap">
  <p class="aperture-cap">The course as a {len(MODS)}-element array, one element per module, grouped into five parts. Pick a learning path below to light up its active aperture; select an element to open that module.</p>
  {aperture}
</div>

<h2 id="paths">Learning paths</h2>
<div class="prose"><p>Each path is an ordered subset of modules for one goal. Earlier modules are prerequisites for later ones, so read them in the order shown.</p></div>
<div class="pathbar" id="pathbar" hidden>{"".join(btns)}</div>
<div class="paths" id="paths-list">{"".join(paths)}</div>

<h2 id="map">The five parts</h2>
<div class="prose"><p>{len(MODS)} modules, {n_demos} interactive demos, {n_eq} numbered equations and {n_scripts} MATLAB/Octave scripts. Each module follows the same arc: derivations with live figures, a visual summary, how experts think about it, a cheat sheet, interview questions, a MATLAB lab, a self-check, figures from the literature and further reading.</p></div>
<div class="map">{"".join(tiers)}</div>

<h2 id="different">What makes it different</h2>
<div class="diff">
  <div><h3>Every figure says where it comes from</h3><p>Each plot carries one of four badges: <span class="prov exact">Exact formula</span><span class="prov sim">Simulation</span><span class="prov model">Approx. model</span><span class="prov schem">Schematic</span>. Simulations name the solver that actually ran; approximate models state their assumptions beside the figure.</p></div>
  <div><h3>Runnable MATLAB and Octave</h3><p>{n_scripts} scripts in base MATLAB with no toolboxes. Every one runs unchanged in GNU Octave 8 and is executed headless on each change; the numbers it prints match the page text.</p></div>
  <div><h3>Literature you can check</h3><p>{n_lit} "Figures from the literature" entries point to a specific figure in a published paper, with citation and DOI, each looked up before it was included. Unverifiable details are removed or listed on the <a href="errata.html">errata page</a>.</p></div>
  <div><h3>Checked, and still learning</h3><p>Every number on these pages was checked by running code. Before release every module also had a careful second review, in which a separate AI model re-derived the key results and looked up the citations, and the claims it could not confirm were rechecked against primary sources. The few still open are listed on the <a href="errata.html">errata page</a>. The most valuable review of all comes from readers, so corrections are always welcome.</p></div>
</div>

<h2 id="sheet">Formula sheet</h2>
<div class="prose"><p>All the key equations of the course, with symbols defined, on one page: <a href="formula-sheet/index.html">open the formula sheet</a>.</p></div>

<h2 id="cite">How to cite</h2>
<div class="prose"><p>Nguyen, M. ({DATE[:4]}). <em>Open Ultrasound Lab</em> (Version {esc(VERSION)}) [Interactive educational resource]. {esc(ARCH)}</p>
<p>BibTeX and CITATION.cff are on the <a href="cite.html">citation page</a>. Found a mistake? <a href="{attr(REPO)}/issues/new">Open an issue</a> or read <a href="contribute.html">how to contribute</a>. Please read the <a href="disclaimer.html">disclaimer</a> before relying on anything here.</p></div>
'''
    static_page('index.html', 'Open Ultrasound Lab — interactive ultrasound signal processing, physics and imaging',
                'A free, open, interactive course: 33 modules on the math, physics, hardware and signal processing of ultrasound imaging, with live figures and runnable MATLAB/Octave.',
                body + js, 'home')


landing()


# ---------------------------------------------------------------- errata
AUDIT_JSON = root / 'audit_summary.json'   # public repo: the audit reports themselves stay private; this summary carries what the site shows


def audit_summary():
    """{'audited': [module numbers], 'flags': [...]} from AUDIT_*.md when present (private source; also refreshes
    audit_summary.json), otherwise from audit_summary.json (public repository)."""
    files = list(root.glob('AUDIT_m*.md'))
    if files:
        flags = audit_flags_md()
        rc = root / 'notes' / 'recheck_2026-10-26.json'   # private: pre-release recheck of the flagged claims
        stats = None
        if rc.exists():
            R = json.loads(rc.read_text()); byc = {r['claim']: r for r in R}
            stats = {k: sum(1 for r in R if r['status'] == k) for k in ('confirmed', 'corrected', 'soften', 'unverifiable', 'open')}
            keep = []
            for fl in flags:
                r = byc.get(fl['claim'])
                if r is None or r['status'] in ('unverifiable', 'open'):
                    if r is not None:
                        fl['verdict'] = 'open'; fl['action'] = r.get('open_note') or 'no accessible source found in the 26 Oct 2026 recheck'
                    keep.append(fl)
            flags = keep
        d = {'audited': sorted({int(m.group(1)) for f in files for m in [re.search(r'm(\d+)', f.stem)] if m}), 'flags': flags, 'recheck': stats}
        AUDIT_JSON.write_text(json.dumps(d, ensure_ascii=False, indent=1))
        return d
    if AUDIT_JSON.exists():
        return json.loads(AUDIT_JSON.read_text())
    return {'audited': [], 'flags': []}


def audit_flags_md():
    """Rows of AUDIT_*.md claim tables whose verdict is unverifiable or disputed."""
    out = []
    for f in sorted(root.glob('AUDIT_*.md'), key=lambda p: [int(x) if x.isdigit() else x for x in re.split(r'(\d+)', p.name)]):
        m = re.search(r'm(\d+)', f.stem)
        n = int(m.group(1)) if m else None
        for line in f.read_text().splitlines():
            if not line.startswith('|') or set(line) <= set('|-: '):
                continue
            cells = [c.strip() for c in re.split(r'(?<!\\)\|', line.strip().strip('|'))]
            vi = next((i for i, c in enumerate(cells) if i >= 1 and re.match(r'(\*\*)?(unverifiable|disputed)\b', c, re.I)), None)
            if vi is None:
                continue
            ci = 1 if cells[0].isdigit() else 0
            if cells[ci].strip().lower() in ('module', 'claim', '#'):
                continue   # table header rows
            verdict = cells[vi].replace('**', '')
            out.append(dict(module=n, claim=cells[ci].replace('\\|', '|'), where=cells[ci + 1] if ci + 1 < vi else '',
                            verdict=verdict, action=cells[-1] if cells[-1] not in ('-', '–', '—') and len(cells) - 1 > vi else ''))
    return out



def glossary_page():
    if not GLOSS:
        return
    items, letters, first = [], [], {}
    for g in sorted(GLOSS, key=lambda g: g['t'].upper()):
        L = g['t'][0].upper()
        if L not in letters:
            letters.append(L); first[L] = gl_id(g['t'])
        href, mn = gl_link(g['t'], 0)
        items.append(f'<li id="{gl_id(g["t"])}" data-k="{attr((g["t"] + " " + g["f"]).lower())}"><span class="t">{esc(g["t"])}</span> '
                     f'<span class="f">{esc(g["f"])}</span><p>{esc(g["d"])}</p>'
                     f'<a href="{attr(href)}">Module {mn:02d} · {esc(MODS[mn]["short"])}</a></li>')
        SEARCH_INDEX['e'].append(['gl', f'{g["t"]}: {g["f"]}. {g["d"]}', mn, GL_FIRST.get(g['t'], (mn, ''))[1]])
    body = (f'<div class="hero"><h1>Glossary</h1><p class="lede">{len(GLOSS)} abbreviations used across the modules. '
            'In the modules, the first use of a term in each section is underlined with dots: hover, tap or press Enter on it to see the definition.</p></div>'
            '<p><input id="gl-filter" type="search" placeholder="Filter terms, e.g. SNR or coherence" aria-label="Filter glossary"></p>'
            '<p class="gl-letters">' + ' '.join(f'<a href="#{first[L]}">{L}</a>' for L in letters) + '</p>'
            '<ul class="gl-list">' + '\n'.join(items) + '</ul>'
            '<script>(function(){var f=document.getElementById("gl-filter"),li=[].slice.call(document.querySelectorAll(".gl-list li"));'
            'f.addEventListener("input",function(){var q=f.value.trim().toLowerCase();li.forEach(function(e){e.hidden=q&&e.dataset.k.indexOf(q)<0;});});})();</script>')
    static_page('glossary.html', 'Glossary — Open Ultrasound Lab', 'Abbreviations and key terms used in Open Ultrasound Lab, with short definitions and links to the modules that teach them.', GL_CSS + body, 'glossary')

def errata_page():
    data = json.loads((root / 'errata.json').read_text())
    flags = AUD['flags']
    audited = AUD['audited']
    ents = []
    for e in data['entries']:
        h = (f'<div class="errata-entry" id="v{attr(e["version"])}"><h3>Version {esc(e["version"])} · {esc(e["date"])}</h3>'
             f'<p><strong>{esc(e["title"])}</strong></p><div class="prose"><p>{e["body"]}</p></div>')
        for it in e.get('items', []):
            h += f'<p>{it}</p>'
        if e.get('flagged') == 'from-audits':
            if audited:
                mods = ', '.join(f'<a href="{MODS[n]["url"]}">{n:02d}</a>' for n in audited if n in MODS)
                h += f'<p>Modules audited: {mods}.</p>'
            if flags:
                lis = []
                for fl in flags:
                    n = fl['module']
                    link = f'<a href="{MODS[n]["url"]}">Module {n:02d}</a>' if n in MODS else 'Module ?'
                    lis.append(f'<li><span class="where">{link}{" · " + esc(fl["where"]) if fl["where"] else ""}</span><br>{esc(fl["claim"])}'
                               f'<br><span class="muted">{esc(fl["verdict"])}{"; " + esc(fl["action"]) if fl["action"] else ""}</span></li>')
                h += (f'<h3 id="flagged">Open checks ({len(flags)})</h3><p class="muted">Claims still awaiting an accessible primary source. '
                      'They may well be correct; treat them with extra care, and if you can confirm or refute one, please tell us.</p>'
                      f'<ol class="flag-list">{"".join(lis)}</ol>')
            else:
                h += ('<h3 id="flagged">Flagged as unverified</h3><ol class="flag-list"><li class="muted">No audit reports are '
                      'present in this build.</li></ol>')
        ents.append(h + '</div>')
    body = (f'<div class="hero"><h1>Errata</h1><p class="lede">Corrections and known issues, newest first. Each entry names the '
            f'version it applies to. To report a new error, use the "Report an error" link on any module page or '
            f'<a href="{attr(REPO)}/issues/new">open an issue</a>.</p></div>' + ''.join(ents))
    static_page('errata.html', 'Errata — Open Ultrasound Lab', 'Corrections, known issues and claims flagged as unverified in Open Ultrasound Lab.', body, 'errata')
    return len(flags)


AUD = audit_summary()
n_flags = errata_page()


# ---------------------------------------------------------------- about, disclaimer, cite, contribute
audited_list = AUD['audited']
aud_txt = (f'{len(audited_list)} of {len(MODS)}' if audited_list else 'none yet')

static_page('about.html', 'About — Open Ultrasound Lab', 'What Open Ultrasound Lab is, how it was built with an AI model, and how it is checked.', f'''
<div class="hero"><h1>About</h1><p class="lede">Open Ultrasound Lab is a free course on the mathematics, physics, hardware and signal processing of ultrasound imaging, written for engineers and researchers who want to read, and then extend, current work in localization microscopy, sound-speed estimation and aberration correction.</p></div>
<div class="prose">
<h2>The project</h2>
<p>The course has {len(MODS)} modules in five parts, from complex numbers and the FFT to the three capstones. Each module derives its results, puts {n_demos // len(MODS)}-odd interactive demos next to the equations, collects the rules of thumb experts quote, and ends with interview questions, a MATLAB lab and a self-check. Modules link to each other by number, so you can start anywhere and follow the prerequisites back.</p>
<h2>How it was built</h2>
<p>M. Nguyen, an ultrasound engineer, planned the course, set the scope and the accuracy rules, and reviewed parts of it. Most of the text, demos and scripts were drafted by an AI model (Claude, by Anthropic) working to written accuracy and style rules (summarized in CONTRIBUTING.md in the repository). Because a language model can state wrong things fluently, the build relies on checks rather than trust:</p>
<ul>
<li><strong>Provenance badges.</strong> Every figure is labeled <span class="prov exact">Exact formula</span>, <span class="prov sim">Simulation</span>, <span class="prov model">Approx. model</span> or <span class="prov schem">Schematic</span>. A simulation badge is allowed only for solvers that actually run in the page or in Octave, and approximate models must state their assumptions.</li>
<li><strong>Render tests.</strong> A headless browser opens every page at desktop and phone width in light and dark themes and fails the build on JavaScript errors, horizontal overflow, canvases that stay blank, MathJax errors or unrendered TeX.</li>
<li><strong>Octave tests.</strong> All {n_scripts} MATLAB scripts run headless in GNU Octave; a script that errors fails the build. The printed numbers were compared with the page text when each script was written.</li>
<li><strong>Independent audit.</strong> A second model instance, without the author's context, re-derives each module's key results blind, extracts every quantitative claim and citation into a table, and marks each correct, wrong, imprecise, unverifiable or disputed. Wrong and imprecise items are fixed and the fixes re-verified; every module ({aud_txt}) was audited before release. Claims the audit could not confirm were then rechecked against primary sources; the few still open are listed on the <a href="errata.html">errata page</a>.</li>
<li><strong>Verified literature.</strong> References and figure numbers in "Figures from the literature" were looked up before inclusion; ones that could not be confirmed were dropped.</li>
</ul>
<p>None of this makes the course error-free. It has been only partly reviewed by a human expert. Read the <a href="disclaimer.html">disclaimer</a>, and please report what you find.</p>
<h2>Who made it</h2>
<p>M. Nguyen, with Claude (Anthropic). More about the author on <a href="{attr(AUTHOR)}">LinkedIn</a>. This is a personal project and does not represent any employer.</p>
<h2 id="license">License</h2>
<p>Text and figures: <a href="https://creativecommons.org/licenses/by/4.0/" rel="license">Creative Commons Attribution 4.0 International</a>. Code (page scripts, MATLAB/Octave files, build tools): MIT. See LICENSE and LICENSE-CONTENT in the <a href="{attr(REPO)}">repository</a>.</p>
</div>''', 'about')

static_page('disclaimer.html', 'Disclaimer — Open Ultrasound Lab', 'Open Ultrasound Lab is an educational resource, partly reviewed, developed with an AI model; not for clinical, regulatory or safety decisions.', f'''
<div class="hero"><h1>Disclaimer</h1></div>
<div class="prose">
<p>This is a free educational resource developed by M. Nguyen with substantial help from an AI model (Claude, by Anthropic) and only partly reviewed. It will contain errors. Every figure is labeled as exact formula, in-browser simulation, approximate model, or schematic; none is measured data unless a cited source says so. Not for clinical, regulatory or safety decisions; not a substitute for IEC/FDA/AIUM standards. Personal project; does not represent the views, products or confidential information of any employer or collaborator. Trademarks belong to their owners.</p>
<p class="muted">Claims still awaiting an accessible source are listed on the <a href="errata.html">errata page</a>. The material is provided "as is", without warranty, under the licenses on the <a href="about.html#license">about page</a>.</p>
</div>''')

cff_src = (root / 'CITATION.cff').read_text()
cff_pub = cff_src.replace('GH_USER', CFG.get('GH_USER', 'GH_USER')).replace('10.5281/zenodo.0000000', CFG['DOI'])
_doi_lines = (f"  publisher    = {{Zenodo}},\n  doi          = {{{CFG['DOI']}}},\n" if HAS_DOI else '')
bib = (f"@misc{{nguyen{DATE[:4]}openultrasoundlab,\n"
       f"  author       = {{Nguyen, M.}},\n"
       f"  title        = {{Open Ultrasound Lab: an interactive course on ultrasound signal processing, physics and imaging}},\n"
       f"  year         = {{{DATE[:4]}}},\n"
       f"  version      = {{{VERSION}}},\n" + _doi_lines +
       f"  url          = {{{CFG['COURSE_URL']}}},\n"
       f"  note         = {{Developed with Claude (Anthropic). Text CC BY 4.0, code MIT}}\n}}")
static_page('cite.html', 'Cite — Open Ultrasound Lab', 'How to cite Open Ultrasound Lab: APA, BibTeX and CITATION.cff.', f'''
<div class="hero"><h1>How to cite</h1><p class="lede">If the course helped your work or teaching, cite it as below{'' if HAS_DOI else '. An archived version with a DOI will be added here'}.</p></div>
<div class="prose"><p>Author: M. Nguyen, ORCID <a href="https://orcid.org/{CFG.get("ORCID","")}">{CFG.get("ORCID","")}</a></p></div>
<h2>APA</h2>
<pre class="cite">Nguyen, M. ({DATE[:4]}). Open Ultrasound Lab: An interactive course on ultrasound signal
    processing, physics and imaging (Version {esc(VERSION)}) [Interactive educational resource].
    {esc(ARCH)}</pre>
<h2>BibTeX</h2>
<pre class="cite">{esc(bib)}</pre>
<h2>CITATION.cff</h2>
<p class="prose">The repository carries this file, so GitHub shows a "Cite this repository" button. <a href="CITATION.cff" download>Download CITATION.cff</a>.</p>
<pre class="cite">{esc(cff_pub)}</pre>
<h2>Citing a single module</h2>
<pre class="cite">Nguyen, M. ({DATE[:4]}). Module 31: ULM end to end. In Open Ultrasound Lab (Version {esc(VERSION)}).
    {esc(ARCH)}</pre>
<p class="prose">Figures you reuse are CC BY 4.0: credit "Open Ultrasound Lab, M. Nguyen" with a link. Figures redrawn from published papers also need the original citation given beside them.</p>''')

static_page('contribute.html', 'Contribute — Open Ultrasound Lab', 'How to report errors and contribute to Open Ultrasound Lab.', f'''
<div class="hero"><h1>Contribute</h1><p class="lede">The most useful contribution is a precise error report: which module, which equation or sentence, what is wrong, and a source or a derivation that shows it.</p></div>
<div class="prose">
<h2>Report an error</h2>
<p>Every module page has a "Report an error" link that opens a GitHub issue titled with the module number. Quote the equation number, for example (15.4), or the section number, and say whether the problem is a wrong number, a sign or factor-of-two, a misattributed citation, a figure that misleads, or a script that fails. Confirmed errors are fixed and listed on the <a href="errata.html">errata page</a>.</p>
<h2>Suggest or write content</h2>
<p>Open an issue first for anything larger than a fix. New or changed module content follows the rules in CONTRIBUTING.md: structure, accuracy rules, provenance badges, MathJax markup and the MATLAB rules (base MATLAB only, runs unchanged in Octave 8).</p>
<h2>Quality gates</h2>
<p>Every pull request runs two checks: the render check opens every page in headless Chromium at 1280 px light and 400 px dark and fails on JavaScript errors, overflow, blank canvases, MathJax errors or raw TeX; the Octave job runs every MATLAB script with a time limit. Both must pass. Details are in CONTRIBUTING.md.</p>
<h2>Where</h2>
<p><a href="{attr(REPO)}">Repository</a> · <a href="{attr(REPO)}/issues">Issues</a> · <a href="{attr(REPO)}/pulls">Pull requests</a></p>
</div>''')


# ---------------------------------------------------------------- formula sheet
def formula_sheet():
    cmd = [sys.executable, str(root / 'build_sheet.py')] + (['--local'] if LOCAL else [])
    subprocess.run(cmd, check=True, capture_output=True)
    s = (root / 'build' / 'formula-sheet.html').read_text()
    if LOCAL:  # serve the shared copy instead of the embedded 2 MB script
        s = re.sub(r'<script>[^<]{100000,}?</script>', f'<script async src="../vendor/{MJ_FILE}"></script>', s, count=1, flags=re.S)

    def kmod(m):
        n = int(m.group(1))
        return f'<h3 class="kmod" id="k-m{n}"><a href="../{MODS[n]["url"]}">' if n in MODS else m.group(0)
    s = re.sub(r'<h3 class="kmod" id="k-m(\d+)"><a href="[^"]*"(?: target="_blank")?(?: rel="noopener")?>', kmod, s)
    s = s.replace('<title>Ultrasound Formula Sheet</title>',
                  '<title>Formula sheet — Open Ultrasound Lab</title>\n<meta name="description" content="Key equations of all 33 Open Ultrasound Lab modules on one page, organized by part.">', 1)
    s = s.replace('</head>', SITE_CSS + '\n' + search_script(1) + '</head>', 1)
    s = s.replace('<body>', '<body>\n' + topbar(1, 'sheet') + '<span id="main"></span>', 1)
    s = s.replace('</body>', f'<div class="page" style="padding-top:0">{footer(1)}</div>\n</body>', 1)
    write('formula-sheet/index.html', s)


# ---------------------------------------------------------------- assets
def assets():
    d = PUB / 'matlab'
    d.mkdir(parents=True, exist_ok=True)
    for f in sorted((root / 'matlab').glob('*.m')):
        shutil.copy2(f, d / f.name)
    (PUB / '.nojekyll').write_text('')
    (PUB / 'CITATION.cff').write_text(cff_pub)
    if not LOCAL:
        shutil.rmtree(PUB / 'vendor', ignore_errors=True)
    if LOCAL:
        (PUB / 'vendor').mkdir(exist_ok=True)
        shutil.copy2(root / 'vendor' / MJ_FILE, PUB / 'vendor' / MJ_FILE)


if __name__ == '__main__':
    pass
formula_sheet()
assets()
glossary_page()
idx_bytes, idx_n = write_search_index()
total = sum(f.stat().st_size for f in PUB.rglob('*') if f.is_file())
print(f'public/: {len(written)} pages, {len(MODS)} modules, {n_questions} interview questions, {n_flags} flagged audit items, '
      f'{link_count} module links rewritten, {total / 1e6:.1f} MB, MathJax {"local" if LOCAL else "CDN"}, search index {idx_n} entries / {idx_bytes / 1e3:.0f} kB')
