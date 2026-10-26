#!/usr/bin/env python3
"""Apply the equal-tuning content to m25.html from tune_out/gen.json."""
import json, pathlib, re
W = pathlib.Path(__file__).parent; P = pathlib.Path('./m25.html')
g = json.loads((W / 'tune_out' / 'gen.json').read_text()); s = P.read_text()
S = g['summary']; LIM = g['lim']; DASF = g['dasF']
def rep(old, new, count=1):
    global s
    n = s.count(old); assert n == count, (n, old[:120]); s = s.replace(old, new)
f3 = lambda v: f'{v:.3f}'
ORDER = ['das', 'mv', 'esbmv', 'cf', 'gcf', 'pcf', 'wiener', 'dmas', 'slsc', 'fxpf', 'admire']

# 1. placeholders
for k in ORDER:
    rep(f'<span id="m25-tuned-{k}" class="m25-tuned">see Section on tuning</span>', f'<span id="m25-tuned-{k}" class="m25-tuned">{g["tuned_txt"][k]} <a href="#m25-tuning">(protocol)</a></span>')
    rep(f'<span id="m25-tuned-tbl-{k}" class="m25-tuned">see Section on tuning</span>', f'<span id="m25-tuned-tbl-{k}" class="m25-tuned">{g["tbl_txt"][k]}</span>')

# 2. scope callout
rep('Each method runs with published or stated defaults (cards above), and will also be shown tuned by one common rule: maximize gCNR subject to point FWHM ≤ 1.1 × DAS.',
    'Each method runs with published or stated defaults (cards above), and is also shown tuned by one common rule: maximize gCNR subject to point FWHM ≤ 1.1 × DAS (the "Settings" selector in Demo 25.6; protocol and trade-off curves in <a href="#m25-tuning">25.7</a>).')

# 3. parameters paragraph + metrics definition of n/a
rep('4-sample windows, 3 bins (25.5, table there).</p>',
    '4-sample windows, 3 bins (25.5, table there). These are the defaults; the "Settings" selector switches to the values tuned by one common rule (<a href="#m25-tuning">25.7</a>; tuned column of the table below).</p>')
rep('pair dip = depth of the valley between the 0.5 mm pair relative to the weaker peak (0 dB = not resolved);',
    'pair dip = depth of the valley between the 0.5 mm pair relative to the weaker peak (0 dB = not resolved; "n/a" = no resolvable pair: the weaker peak or the valley is exactly 0, or the two equal-amplitude peaks differ by more than 12 dB, which happens when noise dominates);')

# 4. Demo 25.6 controls
rep('''        <div class="btn-row"><button class="btn solid" id="m25-g-run" type="button">Run with these settings</button></div>''',
    '''        <div class="ctrl"><label for="m25-g-set">Settings</label><select id="m25-g-set"><option value="default" selected>Default (published or stated values)</option><option value="tuned">Tuned by the common rule (25.7)</option></select></div>
        <div class="readout" id="m25-g-set-o">default settings</div>
        <div class="btn-row"><button class="btn solid" id="m25-g-run" type="button">Run with these settings</button></div>''')
rep('Run reduced ADMIRE live (~45 s)</button>', 'Run reduced ADMIRE live (~@@LIVE_DEF@@ s)</button>')
rep('In this mode the FWHM and pair-dip columns of the table measure noise, not resolution:',
    '<strong>Settings:</strong> "Tuned" reruns every method with the setting chosen by the common rule of <a href="#m25-tuning">25.7</a> (MV, ESBMV, CF, PCF, SLSC and FXPF change; DAS, GCF, Wiener and F-DMAS keep their defaults); the run takes about as long, the 25.7b/c sweeps and the precomputed full-model ADMIRE panel stay at the defaults, and "Run reduced ADMIRE live" then uses the tuned c<sub>λ</sub>. Measured in headless Chromium (our test machine), the live reduced ADMIRE run takes @@LIVE_DEF@@ s with the default c<sub>λ</sub> and @@LIVE_TUN@@ s with the tuned c<sub>λ</sub> (a smaller λ leaves more nonzero coefficients, so the solver needs more sweeps). In line-by-line mode the FWHM and pair-dip columns of the table measure noise, not resolution:')

# 5. ADMIRE card cost text
rep('live reduced model ≈ 45 s per image.', 'live reduced model @@LIVE_DEF@@ s per image (headless Chromium, our test).')

# 6. tuning section after Demo 25.7b
SM = S
txt_sec = f'''  <div class="prose" id="m25-tuning">
    <h3 class="m25-h3">Equal tuning: one rule for every method</h3>
    <p>Defaults are open to the objection that a method was left untuned. So every method is also tuned by one rule, the same for all: <strong>maximize gCNR subject to point FWHM ≤ 1.1 × DAS FWHM</strong>, at −20 dB channel noise and no phase screen. gCNR is the objective because the display cannot move it (Demo 25.7a); the FWHM ceiling stops a method from buying contrast with resolution. The "Settings" selector of Demo 25.6 runs either set, and the methods cards and the summary table give the tuned values.</p>
  </div>
  <div class="callout m25-scope" id="m25-protocol"><span class="label">Tuning protocol</span>
    <p><strong>Objective:</strong> gCNR of the linear output, cyst ROI against background ROI (definitions of 25.6). <strong>Constraint:</strong> lateral −6 dB FWHM of the point at (−3, 24) mm ≤ 1.1 × the FWHM of default DAS (rectangular, 64 elements) on the same data, 1.1 × {DASF:.3f} = {LIM:.3f} mm. Settings whose half-maximum crossings fall outside the ±2.5 mm profile (SLSC with M = 1, three settings) have no valid FWHM and are excluded.</p>
    <p><strong>Data:</strong> a held-out realization of the same phantom, with new scatterer positions and amplitudes (seed 1007) and new channel noise (seed 1108); the displayed run uses seeds 7 and 108. ROI and profile pixels only, as in Demo 25.7b. The chosen setting is then evaluated once on the displayed data (full image, the metrics of the Demo 25.6 table). Everything ran offline in Node with this page's own engine code.</p>
    <p><strong>Grids:</strong> coarse, {g["ncand"]} settings in total, listed in the table below (MV L, K and Δ jointly; ESBMV δ with L; CF exponent including "none"; GCF M<sub>0</sub>; PCF γ; SLSC M with the axial kernel; FXPF p with W; DAS rectangular or Hann apodization). Wiener and F-DMAS have no free parameter in our implementation; FXPF μ, the angles and the aperture stay fixed. <strong>Selection:</strong> the feasible setting with the highest gCNR (exact ties: the narrower point); the default is kept if it is feasible and within 0.005 of that maximum. For scale, the same DAS gives gCNR {f3(SM["das"]["def_ho"])} on the held-out data and {f3(SM["das"]["disp_def"])} on the displayed data: realization-to-realization differences are much larger than 0.005, so gains of a few thousandths below are within that spread.</p>
    <p><strong>ADMIRE:</strong> only c<sub>λ</sub> ∈ {{0.0189/4, 0.0189/2, 0.0189, 0.0189·2}} is searched, with α fixed at the published 0.9, and <strong>only with the reduced model</strong> (≈ 900 predictors per fit). Tuning the full model (≈ 31 000 predictors per fit; its 24-setting grid in Demo 25.5c alone took 3.8 CPU-hours) over a parameter grid was too costly here, so the precomputed full-model panel keeps c<sub>λ</sub> = 0.0189.</p>
    <p><strong>Other objectives would change the winner.</strong> Maximizing CR rewards dynamic-range stretch (ESBMV's default reaches 70.5 dB on the displayed data), minimizing FWHM rewards MV's SNR-dependent point width (25.2), and a different ceiling, noise level or phase screen moves every choice. The trade-off curves below show each method's whole range so that another rule can be applied by eye.</p>
  </div>
'''
FIG = json.dumps(g['fig'], ensure_ascii=False, separators=(',', ':'))
das_h = next(p for p in g['fig']['m']['das']['p'] if 'Hann' in p[2])
pcf = {p[2]: p for p in g['fig']['m']['pcf']['p']}
cf = g['fig']['m']['cf']['p']
demo = f'''  <div class="demo">
    <div class="demo-head"><span class="tag">Demo 25.7d</span><h3>Trade-off curves: point FWHM against gCNR for every setting tried</h3><span class="prov sim">Simulation</span><span class="prov-note">Computed offline in Node with this page's own engine (the M16E and ADM code of this page): {g["ncand"]} settings of 11 methods on the held-out realization (scatterer seed 1007, noise seed 1108), −20 dB channel noise, no phase screen, ROI and profile pixels as in Demo 25.7b; ADMIRE with the reduced model. The values are embedded in the page, not recomputed in your browser.</span></div>
    <div class="demo-body">
      <div class="stage">
        <canvas id="m25-t-cv" aria-label="Point FWHM versus gCNR for every setting of every method, with Pareto fronts, the FWHM ceiling and default and tuned markers"></canvas>
        <div class="legend" id="m25-t-leg"></div>
      </div>
      <div class="controls">
        <div class="ctrl"><label for="m25-t-set">Methods</label><select id="m25-t-set"><option value="all" selected>All eleven</option><option value="adaptive">DAS, MV, ESBMV</option><option value="weights">DAS and the coherence weights</option><option value="other">DAS, F-DMAS, SLSC, FXPF, ADMIRE</option></select></div>
        <div class="ctrl"><label for="m25-t-y">gCNR axis</label><select id="m25-t-y"><option value="zoom" selected>0.85 to 1 (lower settings as ▼ on the bottom edge)</option><option value="full">0 to 1</option></select></div>
        <div class="readout" id="m25-t-read"></div>
      </div>
    </div>
    <div class="demo-foot"><strong>What it shows:</strong> every setting tried, as point FWHM (log axis) against gCNR on the held-out data. Up and to the left is better on both. The line through each method's dots is its Pareto front (no other setting of that method is both narrower and higher), the dashed vertical line is the ceiling 1.1 × DAS FWHM = {LIM:.3f} mm, the open marker is the default and the filled marker the tuned setting (both coincide where the default was kept). Read off: DAS with Hann apodization gains gCNR ({das_h[1]:.3f}) only by widening the point to {das_h[0]:.2f} mm; CF's gCNR falls steadily as the exponent grows ({cf[0][1]:.3f} without the weight, {cf[-1][1]:.3f} at exponent 2, while the point narrows from {cf[0][0]:.3f} to {cf[-1][0]:.3f} mm), so under this rule its best setting is no weight; PCF collapses above γ = 1 ({pcf["γ = 1.25"][1]:.3f}, {pcf["γ = 1.5"][1]:.3f}, {pcf["γ = 2"][1]:.3f} at γ = 1.25, 1.5, 2); SLSC's best gCNR settings all have points wider than 0.6 mm, and only M = 32 with the ±0.30 mm kernel meets the ceiling. <strong>The ranking depends on the objective:</strong> without the ceiling SLSC and ESBMV tie at the top of the held-out data ({f3(max(p[1] for p in g["fig"]["m"]["slsc"]["p"]))}); with it ESBMV leads; a rule on CR, on the narrowest point or on a fixed display range would reorder the methods again.</div>
  </div>
  <div class="tbl-wrap"><table>
    <thead><tr><th>Method</th><th>Grid searched (settings)</th><th class="num">Default: FWHM (mm) / gCNR, held-out</th><th>Tuned setting</th><th class="num">Tuned: FWHM (mm) / gCNR, held-out</th><th class="num">gCNR on the displayed data, default → tuned</th><th>Ceiling binding?</th></tr></thead>
    <tbody>{"".join(g["rows"])}</tbody></table></div>
  <p class="m25-note muted">"Binding" = the highest-gCNR setting of the grid fails the FWHM ceiling, so the constraint decided the choice. ✗ = the default fails the ceiling. ADMIRE rows use the reduced model on both data sets; the gCNR of the precomputed full model on the displayed data is in the Demo 25.6 table. The tuned c<sub>λ</sub> of ADMIRE is better on the held-out data and worse on the displayed data ({f3(SM["admire"]["disp_def"])} → {f3(SM["admire"]["disp_tun"])}), and SLSC must give up gCNR ({f3(SM["slsc"]["disp_def"])} → {f3(SM["slsc"]["disp_tun"])}) to meet the ceiling: the rule is a constraint, not a promise of improvement.</p>
'''
anchor = '''  <div class="prose">
    <p><strong>Aberration robustness.</strong>'''
rep(anchor, txt_sec + demo + anchor)

# 7. JS: tuned constants, groups, runCols, runAll, labels, pair dip, ADMIRE live
rep('''  const parOf = o => Object.assign({}, E.DEF, o || {});
  function newCube(nz)''', f'''  const parOf = o => Object.assign({{}}, E.DEF, o || {{}});
  /* Tuned settings: output of the equal-tuning protocol of 25.7 ("Tuning protocol" box, id m25-protocol). Computed offline in Node with this page's own
     engine (M16E, ADM) on a held-out realization (scatterer seed 1007, noise seed 1108; the displayed run uses 7 and 108) at -20 dB channel noise and no
     phase screen: maximize gCNR subject to point FWHM <= 1.1 x DAS FWHM; default kept if feasible and within 0.005 of the best. {{}} = default kept.
     TUNED_CLAM: ADMIRE c_lambda tuned with the reduced model only (used by the live reduced run in tuned mode). */
  const TUNED = {json.dumps(g["TUNED"])};
  const TUNED_CLAM = {g["TUNED_CLAM"]};
  const ADM_LIVE_S = {{ default: @@LIVE_DEF@@, tuned: @@LIVE_TUN@@ }};   // measured run time (s) of the live reduced ADMIRE run, headless Chromium, our test machine
  const isTunedK = (k, mode) => mode === 'tuned' && !!TUNED[k] && Object.keys(TUNED[k]).length > 0;
  function parGroups(ids, mode) {{   // default: one parameter set; tuned: one group per distinct parameter set
    if (mode !== 'tuned') return parOf(); const gm = new Map();
    ids.forEach(k => {{ const p = parOf(TUNED[k]), key = JSON.stringify(p); if (!gm.has(key)) gm.set(key, [p, []]); gm.get(key)[1].push(k); }});
    return [...gm.values()];
  }}
  function newCube(nz)''')
rep('''      const r = E.runColumn(cube, par, ids, T, null, rows); ids.forEach(k => out[k].push(r[k]));''',
    '''      groups.forEach(([p, gi]) => { const Tg = {}, r = E.runColumn(cube, p, gi, Tg, null, rows); gi.forEach(k => { out[k].push(r[k]); T[k] = (T[k] || 0) + (Tg[k] || 0); }); });''')
rep('''    let cube = newCube(zf.length); const out = {}, T = {}; ids.forEach(k => out[k] = []); let tD = 0, tN = 0, last = performance.now();''',
    '''    let cube = newCube(zf.length); const out = {}, T = {}, groups = Array.isArray(par) ? par : [[par, ids]]; ids.forEach(k => out[k] = []); let tD = 0, tN = 0, last = performance.now();''')
rep('''ST.lbl = $('m25-g-lbl').checked; ST.snrL = null;''', '''ST.lbl = $('m25-g-lbl').checked; ST.set = $('m25-g-set').value; ST.snrL = null;''')
rep('''const r = await runCols(GRID.xs, GRID.zf, ST.scr, GRID.all, parOf(), null, tok,''', '''const PG = parGroups(GRID.all, ST.set), r = await runCols(GRID.xs, GRID.zf, ST.scr, GRID.all, PG, null, tok,''')
rep('''const pa = await runCols(PROF.xs, PROF.zf, ST.scr, GRID.all, parOf(), [16, 16], tok, null, null, L);''', '''const pa = await runCols(PROF.xs, PROF.zf, ST.scr, GRID.all, PG, [16, 16], tok, null, null, L);''')
rep('''const pb = await runCols(PROF.qx, PROF.zf, ST.scr, GRID.all, parOf(), [16, 16], tok, null, null, L);''', '''const pb = await runCols(PROF.qx, PROF.zf, ST.scr, GRID.all, PG, [16, 16], tok, null, null, L);''')
# panel labels
rep('''ctx.fillText(ST.lbl ? STY(c)[k][2] + (lin ? ' (linear)' : '') + ' · line-by-line noise' : isTex ? TEX[mode] + (ST.frame ? ' #' + (ST.frame + 1) : '') : STY(c)[k][2] + (lin ? ' (linear)' : ''), box.x + 6, box.y + 15);''',
    '''const tg = isTunedK(k, ST.set) ? ' · tuned' : ''; ctx.fillText(ST.lbl ? STY(c)[k][2] + tg + (lin ? ' (linear)' : '') + ' · line-by-line noise' : isTex ? TEX[mode] + (ST.frame ? ' #' + (ST.frame + 1) : '') : STY(c)[k][2] + tg + (lin ? ' (linear)' : ''), box.x + 6, box.y + 15);''')
# pair dip
rep('''  function pairDip(qx, prof) { let p1 = 0, p2 = 0, mn = Infinity; qx.forEach((x, i) => { if (Math.abs(x - 2.25) < 0.1) p1 = Math.max(p1, prof[i]); if (Math.abs(x - 2.75) < 0.1) p2 = Math.max(p2, prof[i]); if (x > 2.3 && x < 2.7) mn = Math.min(mn, prof[i]); }); return 20 * Math.log10(Math.max(mn, 1e-12) / Math.max(1e-12, Math.min(p1, p2))); }''',
    '''  /* pair dip (dB, valley re weaker peak); null = no resolvable pair: weaker peak or valley exactly 0 (log undefined), or the two equal-amplitude
     peaks differ by more than 12 dB (noise-dominated). Same rule for every method, including the precomputed ADMIRE row (ADM_DIP_NA). */
  function pairDip(qx, prof) { let p1 = 0, p2 = 0, mn = Infinity; qx.forEach((x, i) => { if (Math.abs(x - 2.25) < 0.1) p1 = Math.max(p1, prof[i]); if (Math.abs(x - 2.75) < 0.1) p2 = Math.max(p2, prof[i]); if (x > 2.3 && x < 2.7) mn = Math.min(mn, prof[i]); });
    const pk = Math.min(p1, p2); if (!(pk > 0) || !(mn > 0) || pk < 0.25 * Math.max(p1, p2)) return null; return 20 * Math.log10(mn / pk); }
  /* precomputed full-model ADMIRE settings (noise dB | screen ns) where pairDip() returns null, evaluated offline on the stored pair profiles with the same rule */
  const ADM_DIP_NA = ['30|20', '20|40', '30|40'];
  const dipCell = (v, f) => v === null ? '<span title="no resolvable pair">n/a</span>' : f(v, 1);''')
rep('''const f = (v, d) => v === undefined ? '…' : fmt(v, d);
      h += '<tr><td><span style="display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;background:' + s[k][0] + '"></span>' + s[k][2] + '</td><td class="num">' + f(m.fwhm, 3) + '</td><td class="num">' + (m.psl !== undefined && m.psl < -150 ? 'none' : f(m.psl, 1)) + '</td><td class="num">' + f(m.dip, 1) + '</td>''',
    '''const f = (v, d) => v === undefined ? '…' : fmt(v, d);
      h += '<tr><td><span style="display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;background:' + s[k][0] + '"></span>' + s[k][2] + (isTunedK(k, ST.set) ? ' <span class="muted">(tuned)</span>' : '') + '</td><td class="num">' + f(m.fwhm, 3) + '</td><td class="num">' + (m.psl !== undefined && m.psl < -150 ? 'none' : f(m.psl, 1)) + '</td><td class="num">' + dipCell(m.dip, f) + '</td>''')
rep('''    GRID.all.forEach(k => { const el = $('m25-mc-t-' + k), t = ST.times && ST.times[k]; if (el && t !== undefined) el.textContent''',
    '''    drawSetO(); if (ST.set !== 'tuned') GRID.all.forEach(k => { const el = $('m25-mc-t-' + k), t = ST.times && ST.times[k]; if (el && t !== undefined) el.textContent''')
rep('''  function drawTable() {''', '''  function drawSetO() {   // settings readout: what the shown run used, and whether the selector differs
    const o = $('m25-g-set-o'); if (!o) return; const want = $('m25-g-set').value, used = ST.img ? ST.set : null;
    const ch = GRID.all.filter(k => isTunedK(k, 'tuned')).map(k => STY(Lab.colors())[k][2]).join(', ');
    const d = m => m === 'tuned' ? 'tuned settings (' + ch + ' changed; the others keep their defaults; ADMIRE live uses c_λ = 0.0189/4)' : 'default settings';
    o.textContent = (used ? 'shown: ' + d(used) : d(want)) + (used && used !== want ? '; press Run to switch to ' + (want === 'tuned' ? 'tuned' : 'default') + ' settings' : '');
  }
  function drawTable() {''')
rep('''  $('m25-g-run').addEventListener('click', () => runAll());''', '''  $('m25-g-run').addEventListener('click', () => runAll());
  $('m25-g-set').addEventListener('change', () => { drawSetO(); drawAdmPanel(); });''')
# ADMIRE live with tuned c_lambda
rep('''    const tok = ST.token, par = ADM.RED, S = ADM.setup(par),''', '''    const tok = ST.token, clam = ST.set === 'tuned' ? TUNED_CLAM : ADM.RED.clam, par = clam === ADM.RED.clam ? ADM.RED : Object.assign({}, ADM.RED, { clam }), S = ADM.setup(par),''')
rep('''    const L = ST.admLive = { tok, busy: true, cols: [], prog: 'line 0 / ' + GRID.xs.length };''', '''    const L = ST.admLive = { tok, busy: true, cols: [], clam, prog: 'line 0 / ' + GRID.xs.length };''')
rep('''      return { data, t2: 'reduced · live', t3: ST.lbl ? 'shared-recording noise' : '',
        say: 'ADMIRE panel: reduced model (≈ 900 predictors per fit), computed live on this run\\'s data (' ''',
    '''      const tc = L.clam !== ADM.RED.clam;
      return { data, t2: 'reduced · live' + (tc ? ' · tuned c_λ' : ''), t3: ST.lbl ? 'shared-recording noise' : '',
        say: 'ADMIRE panel: reduced model (≈ 900 predictors per fit), c_λ = ' + (tc ? '0.0189/4 (tuned)' : '0.0189 (default)') + ', computed live on this run\\'s data (' ''')
rep('''        + '. "Run reduced ADMIRE live" computes the reduced model for the exact settings.' };''',
    '''        + '. "Run reduced ADMIRE live" computes the reduced model for the exact settings' + (ST.set === 'tuned' ? ' with the tuned c_λ = 0.0189/4 (this full-model panel keeps the default 0.0189).' : '.') };''')
rep('''    const btn = $('m25-g-adm-run'); if (btn) btn.disabled = !ST.sim || !ST.img || ST.busy || !!(ST.admLive && ST.admLive.busy);''',
    '''    const btn = $('m25-g-adm-run'); if (btn) { btn.disabled = !ST.sim || !ST.img || ST.busy || !!(ST.admLive && ST.admLive.busy); const tm = ST.img && ST.set === 'tuned' ? 'tuned' : 'default'; btn.textContent = 'Run reduced ADMIRE live (~' + ADM_LIVE_S[tm] + ' s' + (tm === 'tuned' ? ', tuned c_λ' : '') + ')'; }''')
rep('''      return { lbl: 'live, reduced', cr: rm.cr,''', '''      return { lbl: 'live, reduced' + (L.clam !== ADM.RED.clam ? ', tuned c_λ' : ''), cr: rm.cr,''')
rep('''fwhm: m.fwhm, psl: m.psl, dip: m.dip, zc: m.zc,''', '''fwhm: m.fwhm, psl: m.psl, dip: ADM_DIP_NA.includes(g.nd + '|' + g.ab) ? null : m.dip, zc: m.zc,''')
rep('''+ f(m.fwhm, 3) + '</td><td class="num">' + (m.psl !== undefined && m.psl !== null && m.psl < -150 ? 'none' : f(m.psl, 1)) + '</td><td class="num">' + f(m.dip, 1) + '</td>''',
    '''+ f(m.fwhm, 3) + '</td><td class="num">' + (m.psl !== undefined && m.psl !== null && m.psl < -150 ? 'none' : f(m.psl, 1)) + '</td><td class="num">' + dipCell(m.dip, f) + '</td>''')

# 8. trade-off figure JS (after the 25.7b listener)
fig_js = '''
  /* ======================= Demo 25.7d: trade-off curves of the equal-tuning protocol (precomputed offline, see TUNED) ======================= */
  /* T7C: every candidate [FWHM (mm), gCNR, setting] on the held-out realization; d / t = index of the default / tuned setting; lim = 1.1 x DAS FWHM */
  const T7C = ''' + FIG + ''';
  const T_ALL = ['das', 'mv', 'esbmv', 'cf', 'gcf', 'pcf', 'wiener', 'dmas', 'slsc', 'fxpf', 'admire'];
  const T_SETS = { all: T_ALL, adaptive: ['das', 'mv', 'esbmv'], weights: ['das', 'cf', 'gcf', 'pcf', 'wiener'], other: ['das', 'dmas', 'slsc', 'fxpf', 'admire'] };
  function drawTune() {
    const cv = $('m25-t-cv'); if (!cv) return; const c = Lab.colors(), s = STY(c), ids = T_SETS[$('m25-t-set').value], zoom = $('m25-t-y').value === 'zoom';
    const col = k => k === 'admire' ? c.good : s[k][0], dash = k => k === 'admire' ? [2, 3] : s[k][1], nm = k => k === 'admire' ? 'ADMIRE (reduced)' : s[k][2];
    const { ctx, w, h } = Lab.canvas(cv, 360), lg = Math.log10, yr = zoom ? [0.85, 1.0] : [0, 1.0];
    const ax = Lab.axes(ctx, { x: 52, y: 22, w: w - 66, h: h - 64 }, [lg(0.012), lg(5)], yr, { xticks: [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5].map(lg), xfmt: v => String(+Math.pow(10, v).toPrecision(1)),
      yticks: zoom ? [0.85, 0.9, 0.95, 1] : [0, 0.25, 0.5, 0.75, 1], yfmt: v => v.toFixed(2), xlabel: 'point FWHM (mm, log axis)', ylabel: 'gCNR (held-out data)', title: 'every setting tried · ○ default · ● tuned' });
    const X = v => ax.X(lg(v)), Y = v => ax.Y(Math.max(yr[0], v));
    ctx.save(); ctx.strokeStyle = c.fg; ctx.setLineDash([5, 4]); ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(X(T7C.lim), ax.box.y); ctx.lineTo(X(T7C.lim), ax.box.y + ax.box.h); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = c.fg; ctx.font = '11px ' + c.mono; ctx.textAlign = 'right'; ctx.textBaseline = 'top'; ctx.fillText('1.1 × DAS FWHM = ' + T7C.lim.toFixed(3) + ' mm', X(T7C.lim) - 4, ax.box.y + 4); ctx.restore();
    const mark = (k, p, filled, r) => { const below = p[1] < yr[0], x = X(p[0]), y = Y(p[1]); ctx.beginPath();
      if (below) { ctx.moveTo(x - r, y - 2 * r); ctx.lineTo(x + r, y - 2 * r); ctx.lineTo(x, y); ctx.closePath(); } else ctx.arc(x, y, r, 0, TAU);
      if (filled) { ctx.fillStyle = col(k); ctx.fill(); ctx.strokeStyle = c.bg; ctx.lineWidth = 1; ctx.stroke(); } else { ctx.strokeStyle = col(k); ctx.lineWidth = 2; ctx.stroke(); } };
    ax.clip();
    ids.forEach(k => { const M = T7C.m[k], pts = M.p; ctx.save(); ctx.globalAlpha = 0.45; pts.forEach(p => mark(k, p, true, 2.2)); ctx.restore();
      const front = []; let best = -1; pts.slice().sort((a, b) => a[0] - b[0] || b[1] - a[1]).forEach(p => { if (p[1] > best) { front.push(p); best = p[1]; } });
      if (front.length > 1) { ctx.save(); ctx.strokeStyle = col(k); ctx.lineWidth = 1.8; if (dash(k)) ctx.setLineDash(dash(k)); ctx.beginPath(); front.forEach((p, i) => { const x = X(p[0]), y = Y(p[1]); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }); ctx.stroke(); ctx.restore(); } });
    ids.forEach(k => { const M = T7C.m[k]; mark(k, M.p[M.d], false, 6.5); mark(k, M.p[M.t], true, 4.2); });
    ax.unclip();
    $('m25-t-leg').innerHTML = ids.map(k => '<span><i style="background:' + col(k) + '"></i>' + nm(k) + '</span>').join('') + '<span>dots: settings · line: Pareto front · dashed: FWHM ceiling</span>';
    $('m25-t-read').innerHTML = ids.map(k => { const M = T7C.m[k], d = M.p[M.d], t = M.p[M.t];
      return '<b>' + nm(k) + '</b>: default ' + d[0].toFixed(3) + ' mm / ' + d[1].toFixed(3) + (M.d === M.t ? ' (kept)' : ' → ' + t[2] + ': ' + t[0].toFixed(3) + ' mm / ' + t[1].toFixed(3)); }).join('<br>');
  }
  $('m25-t-set').addEventListener('change', drawTune); $('m25-t-y').addEventListener('change', drawTune);
'''
rep('''  $('m25-s-par').addEventListener('change', () => { if (ST.sim && !ST.busy) { const tok = ST.token; sensRun(tok).catch(e => { if (e !== 'stale') console.error(e); }); } else drawSens(); });''',
    '''  $('m25-s-par').addEventListener('change', () => { if (ST.sim && !ST.busy) { const tok = ST.token; sensRun(tok).catch(e => { if (e !== 'stale') console.error(e); }); } else drawSens(); });''' + fig_js)
rep('''    drawAP(); drawMV(); drawFX(); drawGrid(); drawDR(); drawSens(); drawAb(); drawRules(); drawAM(); drawHW();
    admLoad()''', '''    drawAP(); drawMV(); drawFX(); drawGrid(); drawDR(); drawSens(); drawAb(); drawTune(); drawRules(); drawAM(); drawHW();
    admLoad()''')
rep('''if (Lab.current === 'm25') { drawAP(); drawMV(); drawFX(); drawGrid(); drawDR(); drawSens(); drawAb(); drawRules();''', '''if (Lab.current === 'm25') { drawAP(); drawMV(); drawFX(); drawGrid(); drawDR(); drawSens(); drawAb(); drawTune(); drawRules();''')
P.write_text(s); print('patched')
