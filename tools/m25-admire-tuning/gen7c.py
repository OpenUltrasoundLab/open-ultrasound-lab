#!/usr/bin/env python3
"""Equal tuning: turn tune_out/{classic,admire}.json into page content for m25.html (placeholders, constants, figure data, results table)."""
import json, pathlib, math
W = pathlib.Path(__file__).parent
C = json.loads((W / 'tune_out' / 'classic.json').read_text()); A = json.loads((W / 'tune_out' / 'admire.json').read_text())
LIM = C['lim']; DASF = C['dasFwhm']
ORDER = ['das', 'mv', 'esbmv', 'cf', 'gcf', 'pcf', 'wiener', 'dmas', 'slsc', 'fxpf', 'admire']
NAME = {'das': 'DAS', 'mv': 'MV', 'esbmv': 'ESBMV', 'cf': 'CF', 'gcf': 'GCF', 'pcf': 'PCF', 'wiener': 'Wiener', 'dmas': 'F-DMAS', 'slsc': 'SLSC', 'fxpf': 'FXPF', 'admire': 'ADMIRE'}
ho = dict(C['ho']); ho['admire'] = A['ho']
pick = dict(C['pick']); pick['admire'] = A['pick']
disp = {}
for d in C['disp']: disp.setdefault(d['k'], {})[d['role']] = d
dd = {d['par']['clam']: d for d in A['disp']}
disp['admire'] = {'default': dd[0.0189]}
if not A['pick']['same']: disp['admire']['tuned'] = dd[A['pick']['par']['clam']]

def sup(e): return {'-4': '10<sup>−4</sup>', '-3': '10<sup>−3</sup>', '-2': '10<sup>−2</sup>', '-1': '10<sup>−1</sup>', '0': '1', '1': '10'}[e]
def D(v): return sup(str(round(math.log10(v))))
def ptxt(k, p):
    if k == 'das': return 'rectangular apodization' if p['apod'] == 'rect' else 'Hann apodization'
    if k == 'mv': return f"L = {p['L']}, K = ±{p['K']}, Δ = {D(p['D'])}"
    if k == 'esbmv': return f"L = {p['L']}, δ = {p['delta']:g} (K, Δ as MV default)"
    if k == 'cf': return 'no CF weight (exponent 0, i.e. DAS)' if p['cfp'] == 0 else f"CF<sup>{p['cfp']:g}</sup>"
    if k == 'gcf': return f"M<sub>0</sub> = {p['M0']}"
    if k == 'pcf': return f"γ = {p['gamma']:g}"
    if k == 'slsc': return f"M = {p['M']}, kernel ±{p['Ks'] * 0.01875:.2f} mm"
    if k == 'fxpf': return f"p = {p['p']}, W = {p['W']} ({p['W'] * 0.075:.1f} mm)"
    if k == 'admire':
        r = {0.0189 / 4: '0.0189/4', 0.0189 / 2: '0.0189/2', 0.0189: '0.0189', 0.0189 * 2: '0.0189·2'}
        return 'c<sub>λ</sub> = ' + next(v for kk, v in r.items() if abs(kk - p['clam']) < 1e-9)
    return 'no free parameter'
def f3(v): return f'{v:.3f}'
def best_unc(k):
    ok = [r for r in ho[k] if r.get('fwhmOK', True)]; g = max(r['gcnr'] for r in ok); return [r for r in ok if r['gcnr'] == g]
def binding(k):
    return all(r['fwhm'] > LIM for r in best_unc(k)) and len(ho[k]) > 1
defrec = {k: next(r for r in ho[k] if r['def']) for k in ORDER}
pkrec = {k: next(r for r in ho[k] if r['par'] == pick[k]['par']) for k in ORDER}

# ---------- placeholder texts ----------
tuned_txt = {}
for k in ORDER:
    d0 = disp[k]['default']
    if k in ('wiener', 'dmas'):
        tuned_txt[k] = 'default already optimal under the rule (no free parameter in our implementation)'; continue
    if pick[k]['same']:
        extra = ''
        if k == 'das':
            h = next(r for r in ho['das'] if r['par']['apod'] == 'hann')
            extra = f'; Hann apodization raises held-out gCNR to {f3(h["gcnr"])} but widens the point to {h["fwhm"]:.2f} mm, beyond the {LIM:.2f} mm limit'
        tuned_txt[k] = 'default already optimal under the rule' + extra; continue
    d1 = disp[k]['tuned']; s = f'tuned: {ptxt(k, pick[k]["par"])} (gCNR {f3(d0["gcnr"])} → {f3(d1["gcnr"])})'
    if k == 'slsc': s += f'; the default point FWHM ({defrec[k]["fwhm"]:.2f} mm) fails the {LIM:.2f} mm limit'
    if k == 'admire':
        s = (f'tuned on the reduced model only: {ptxt(k, pick[k]["par"])} (reduced-model gCNR {f3(d0["gcnr"])} → {f3(d1["gcnr"])} on the displayed data, '
             f'{f3(defrec[k]["gcnr"])} → {f3(pkrec[k]["gcnr"])} on the held-out data; the default point FWHM {defrec[k]["fwhm"]:.4f} mm fails the {LIM:.4f} mm limit)')
    tuned_txt[k] = s
tbl_txt = {}
for k in ORDER:
    if k in ('wiener', 'dmas'): tbl_txt[k] = 'no free parameter'
    elif pick[k]['same']: tbl_txt[k] = 'default already optimal'
    else:
        tbl_txt[k] = f'{ptxt(k, pick[k]["par"])} (gCNR {f3(disp[k]["default"]["gcnr"])} → {f3(disp[k]["tuned"]["gcnr"])}' + (', reduced model' if k == 'admire' else '') + ')'

# ---------- constants for the page ----------
TUNED = {k: (pick[k]['par'] if not pick[k]['same'] else {}) for k in ORDER if k != 'admire'}
TUNED_CLAM = pick['admire']['par']['clam']

# ---------- figure data: every candidate on the held-out data ----------
fig = {'lim': round(LIM, 4), 'dasF': round(DASF, 4), 'm': {}}
for k in ORDER:
    pts = []; di = ti = None; nbad = 0
    for r in ho[k]:
        if not r.get('fwhmOK', True): nbad += 1; continue
        if r['def']: di = len(pts)
        if r['par'] == pick[k]['par']: ti = len(pts)
        pts.append([round(r['fwhm'], 4), round(r['gcnr'], 4), ptxt(k, r['par']).replace('<sup>', '^').replace('</sup>', '').replace('<sub>', '').replace('</sub>', '')])
    fig['m'][k] = {'p': pts, 'd': di, 't': ti, 'bad': nbad}

# ---------- results table ----------
GRIDTXT = {
  'das': 'apodization {rectangular, Hann}',
  'mv': 'L {8, 16, 24, 32, 40, 48} × K {0…4} × Δ {10<sup>−4</sup> … 10} (decades)',
  'esbmv': 'L {16, 24, 32, 40, 48} × δ {0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9}',
  'cf': 'exponent {0 (none), 0.25, 0.5, 1, 1.5, 2}',
  'gcf': 'M<sub>0</sub> {0, 1, 2, 3, 4, 6, 8}',
  'pcf': 'γ {0.25, 0.5, 0.75, 1, 1.25, 1.5, 2}',
  'wiener': 'none (no parameter)', 'dmas': 'none (band-pass fixed at √2 σ<sub>f</sub>)',
  'slsc': 'M {1, 2, 4, 6, 8, 12, 16, 24, 32} × kernel ±{0.075, 0.15, 0.30} mm',
  'fxpf': 'p {1, 2, 3, 4, 6, 8} × W {8, 16, 32} samples',
  'admire': 'c<sub>λ</sub> {0.0189/4, 0.0189/2, 0.0189, 0.0189·2}, reduced model, α = 0.9'}
rows = []
for k in ORDER:
    d, t = defrec[k], pkrec[k]; n = len(ho[k])
    dd0 = disp[k]['default']; dd1 = disp[k].get('tuned', dd0)
    same = pick[k]['same']
    rows.append('<tr><td>' + NAME[k] + '</td><td>' + GRIDTXT[k] + f' <span class="muted">({n})</span></td><td class="num">' + f'{d["fwhm"]:.3f} / {f3(d["gcnr"])}' + ('' if d['fwhm'] <= LIM else ' <span class="muted" title="fails the FWHM limit">✗</span>')
                + '</td><td>' + ('= default' if same else ptxt(k, pick[k]['par'])) + '</td><td class="num">' + ('—' if same else f'{t["fwhm"]:.3f} / {f3(t["gcnr"])}')
                + '</td><td class="num">' + (f3(dd0['gcnr']) if same else f'{f3(dd0["gcnr"])} → {f3(dd1["gcnr"])}') + '</td><td>' + ('yes' if binding(k) else 'no') + '</td></tr>')
summary = {k: {'def_ho': defrec[k]['gcnr'], 'tun_ho': pkrec[k]['gcnr'], 'def_fw': defrec[k]['fwhm'], 'tun_fw': pkrec[k]['fwhm'], 'disp_def': disp[k]['default']['gcnr'], 'disp_tun': disp[k].get('tuned', disp[k]['default'])['gcnr'],
               'disp_def_fw': disp[k]['default']['fwhm'], 'disp_tun_fw': disp[k].get('tuned', disp[k]['default'])['fwhm'], 'binding': binding(k), 'same': pick[k]['same'], 'n': len(ho[k])} for k in ORDER}
out = {'tuned_txt': tuned_txt, 'tbl_txt': tbl_txt, 'TUNED': TUNED, 'TUNED_CLAM': TUNED_CLAM, 'fig': fig, 'rows': rows, 'summary': summary, 'lim': LIM, 'dasF': DASF,
       'ncand': sum(len(ho[k]) for k in ORDER), 'admMs': {str(d['par']['clam']): d['ms'] for d in A['disp']}}
(W / 'tune_out' / 'gen.json').write_text(json.dumps(out, ensure_ascii=False, indent=1))
for k in ORDER: print(k, '|', tuned_txt[k]);
print(json.dumps(summary, indent=0)[:3000]); print('ncand', out['ncand'], 'lim', LIM)
