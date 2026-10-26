// Equal tuning protocol for module 25 (Demo 25.6).
// Rule: maximize gCNR (linear output, cyst vs background ROI) subject to point FWHM <= 1.1 x DAS FWHM (default DAS, rectangular),
// channel noise -20 dB, no phase screen. Tuning data: HELD-OUT realization (scatterer seed HO, noise seed HO+101) != displayed (7 / 108).
// Then every default and tuned setting is evaluated on the displayed realization with the page's own Demo 25.6 metric definitions.
// Uses the page's own engine (M16E + ADM + GRID/PROF/ROI) read from m25.html.
// Usage: node tune.js <part>   part = classic | admire | all
const fs = require('fs'), path = require('path');
const { Lab, loadAll } = require('./load.js'); Lab.nextPow2 = n => 1 << Math.ceil(Math.log2(Math.max(2, n)));   // copied from shell.html
const { E, ADM } = loadAll();
const page = fs.readFileSync('./m25.html', 'utf8');
const gsrc = page.slice(page.indexOf('  const GRID = { ids:'), page.indexOf('ROI.rows = [8, 61];') + 'ROI.rows = [8, 61];'.length);
const { GRID, PROF, ROI } = new Function('E', gsrc + '\nreturn { GRID, PROF, ROI };')(E);
const part = process.argv[2] || 'all', OUT = path.join(__dirname, 'tune_out'); fs.mkdirSync(OUT, { recursive: true });
const HO = 1007, SHOW = 7, NOISE = -20, TOL = 0.005, FWHM_FACTOR = 1.1;
const N = E.NEL, now = () => Number(process.hrtime.bigint()) / 1e6;
const newCube = nz => ({ re: new Float32Array(nz * N), im: new Float32Array(nz * N), nz });
const scr0 = E.screen(0, 3, 33);
function roiSplit(colsArr, xs, z0, dz, rlo, rhi) {   // identical to the page
  const rc = E.scene.roiCyst, rb = E.scene.roiBg, cy = [], bg = [];
  for (let c = 0; c < xs.length; c++) { const x = xs[c], col = colsArr[c];
    for (let r = rlo; r <= rhi; r++) { const z = Math.round((z0 + r * dz) * 1e4) / 1e4, v = Math.max(col[r], 0);
      if ((x - rc[0]) ** 2 + (z - rc[1]) ** 2 < rc[2] ** 2) cy.push(v); if (x >= rb[0] - 1e-9 && x <= rb[1] + 1e-9 && z >= rb[2] - 1e-9 && z <= rb[3] + 1e-9) bg.push(v); } }
  return { cy, bg };
}
// pair dip with a validity test (the version added to the page in 7c): null = no resolvable pair, -Infinity = valley exactly 0
function pairDipM(qx, prof) {
  let p1 = 0, p2 = 0, mn = Infinity, bs = 0, bn = 0;
  qx.forEach((x, i) => { if (Math.abs(x - 2.25) < 0.1) p1 = Math.max(p1, prof[i]); if (Math.abs(x - 2.75) < 0.1) p2 = Math.max(p2, prof[i]); if (x > 2.3 && x < 2.7) mn = Math.min(mn, prof[i]);
    if (Math.abs(x - 2.5) >= 0.6 - 1e-9) { bs += prof[i] * prof[i]; bn++; } });
  const pk = Math.min(p1, p2);
  if (!(pk > 0) || !(mn > 0) || pk < 0.25 * Math.max(p1, p2)) return null; return 20 * Math.log10(mn / pk);
}
function pairDipOld(qx, prof) { let p1 = 0, p2 = 0, mn = Infinity; qx.forEach((x, i) => { if (Math.abs(x - 2.25) < 0.1) p1 = Math.max(p1, prof[i]); if (Math.abs(x - 2.75) < 0.1) p2 = Math.max(p2, prof[i]); if (x > 2.3 && x < 2.7) mn = Math.min(mn, prof[i]); }); return 20 * Math.log10(Math.max(mn, 1e-12) / Math.max(1e-12, Math.min(p1, p2))); }

/* ---------- candidate grids (coarse) ---------- */
const G = {
  das: [{ apod: 'rect' }, { apod: 'hann' }],
  mv: [], esbmv: [],
  cf: [0, 0.25, 0.5, 1, 1.5, 2].map(v => ({ cfp: v })),
  gcf: [0, 1, 2, 3, 4, 6, 8].map(v => ({ M0: v })),
  pcf: [0.25, 0.5, 0.75, 1, 1.25, 1.5, 2].map(v => ({ gamma: v })),
  wiener: [{}], dmas: [{}],
  slsc: [], fxpf: [] };
[8, 16, 24, 32, 40, 48].forEach(L => [0, 1, 2, 3, 4].forEach(K => [1e-4, 1e-3, 1e-2, 1e-1, 1, 10].forEach(D => G.mv.push({ L, K, D }))));
[16, 24, 32, 40, 48].forEach(L => [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9].forEach(delta => G.esbmv.push({ L, delta })));
[1, 2, 4, 6, 8, 12, 16, 24, 32].forEach(M => [4, 8, 16].forEach(Ks => G.slsc.push({ M, Ks })));
[1, 2, 3, 4, 6, 8].forEach(p => [8, 16, 32].forEach(W => G.fxpf.push({ p, W })));
const KEYS = { das: ['apod'], mv: ['L', 'K', 'D'], esbmv: ['L', 'delta'], cf: ['cfp'], gcf: ['M0'], pcf: ['gamma'], wiener: [], dmas: [], slsc: ['M', 'Ks'], fxpf: ['p', 'W'] };
const isDef = (k, c) => KEYS[k].every(q => c[q] === E.DEF[q]);
const ADM_CL = [0.0189 / 4, 0.0189 / 2, 0.0189, 0.0189 * 2];

/* ---------- held-out: ROI + profile cubes (like the page's 25.7b sweeps) ---------- */
function heldOutClassic() {
  const t0 = now(), sim = E.simulate({ seed: HO, ampPoint: 15, ampBright: 60, noiseDB: NOISE, keepNoise: false });
  const roiC = ROI.xs.map(x => { const c = newCube(ROI.zf.length); E.column(sim, x, ROI.zf, scr0, c); return c; });
  const profC = PROF.xs.map(x => { const c = newCube(PROF.zf.length); E.column(sim, x, PROF.zf, scr0, c); return c; });
  console.log('held-out cubes', ((now() - t0) / 1000).toFixed(1), 's');
  const res = {};
  for (const k of Object.keys(G)) {
    res[k] = [];
    for (const cand of G[k]) {
      const par = Object.assign({}, E.DEF, cand), t1 = now(), T = {};
      const a = roiC.map(c => E.runColumn(c, par, [k], T, null, ROI.rows)[k]), b = profC.map(c => E.runColumn(c, par, [k], T, null, [16, 16])[k]);
      const { cy, bg } = roiSplit(a, ROI.xs, ROI.z0, E.DZ, ROI.rows[0], ROI.rows[1]); let mx = 0; a.forEach(col => { for (let r = ROI.rows[0]; r <= ROI.rows[1]; r++) mx = Math.max(mx, col[r]); });
      const prof = b.map(col => Math.max(col[16], 0)), pm = E.profileMetrics(PROF.xs, prof, -3), rm = E.regionMetrics(cy, bg, mx, []);
      res[k].push({ par: cand, def: isDef(k, cand), fwhm: pm.fwhm, fwhmOK: fwhmValid(prof, pm), gcnr: rm.gcnr, cr: rm.cr, ssnr: rm.ssnr, ms: now() - t1 });
    }
    console.log(k, res[k].length, 'candidates', (res[k].reduce((s, r) => s + r.ms, 0) / 1000).toFixed(1), 's');
  }
  return res;
}
function heldOutAdmire() {
  const sim = E.simulate({ seed: HO, ampPoint: 15, ampBright: 60, noiseDB: NOISE, keepNoise: false }), out = [];
  const z0 = GRID.z0, nz = 164, zf = []; for (let i = 0; i <= (nz - 1) * E.SUB; i++) zf.push(z0 + i * E.DZF);
  const cubes = ROI.xs.map(x => { const c = newCube(zf.length); E.column(sim, x, zf, scr0, c); return c; });
  const pz = PROF.zf.slice(0, 19 * E.SUB + 1), pcs = PROF.xs.map(x => { const c = newCube(pz.length); E.column(sim, x, pz, scr0, c); return c; });
  for (const cl of ADM_CL) {
    const par = Object.assign({}, ADM.RED, { clam: cl }), S = ADM.setup(par), t0 = now();
    const cols = cubes.map((c, i) => Array.from(ADM.line(c, ROI.xs[i], z0, S, par).subarray(0, GRID.nz)));
    const prof = pcs.map((c, i) => Math.max(ADM.line(c, PROF.xs[i], PROF.zf[0], S, par, S.hop < S.Lw ? [14, 16] : [14])[16], 0));
    const { cy, bg } = roiSplit(cols, ROI.xs, GRID.z0, E.DZ, 0, GRID.nz - 1); let mx = 0; cols.forEach(col => col.forEach(v => { if (v > mx) mx = v; }));
    const pm = E.profileMetrics(PROF.xs, prof, -3), rm = E.regionMetrics(cy, bg, mx, []);
    out.push({ par: { clam: cl }, def: cl === 0.0189, fwhm: pm.fwhm, fwhmOK: fwhmValid(prof, pm), gcnr: rm.gcnr, cr: rm.cr, ssnr: rm.ssnr, zc: cy.filter(v => v === 0).length / cy.length, ms: now() - t0 });
    console.log('admire clam', cl, JSON.stringify(out[out.length - 1]));
  }
  return out;
}
// feasible: valid FWHM <= lim. Pick = max gCNR (exact ties: smallest FWHM); the default is kept if it is feasible and within TOL of that max.
function choose(list, lim) {
  const d = list.find(r => r.def) || list[0], feas = list.filter(r => r.fwhmOK && r.fwhm <= lim);
  if (!feas.length) { const best = list.filter(r => r.fwhmOK).reduce((b, r) => r.fwhm < b.fwhm ? r : b); return { pick: best, feasible: false, def: d }; }
  const gmax = Math.max(...feas.map(r => r.gcnr));
  if (feas.includes(d) && d.gcnr >= gmax - TOL - 1e-12) return { pick: d, feasible: true, def: d, kept: true };
  return { pick: feas.filter(r => r.gcnr === gmax).reduce((b, r) => r.fwhm < b.fwhm ? r : b), feasible: true, def: d };
}
// FWHM is valid only if both half-maximum crossings lie inside the measured profile (+-2.5 mm around the point) and the width is positive
const fwhmValid = (prof, pm) => isFinite(pm.fwhm) && pm.fwhm > 0 && prof[0] < pm.peak / 2 && prof[prof.length - 1] < pm.peak / 2;

/* ---------- displayed realization: full image + profiles, exactly like Demo 25.6's table ---------- */
function displayedEval(sets) {   // sets: [{k, par}] ; returns metrics per entry
  const sim = E.simulate({ seed: SHOW, ampPoint: 15, ampBright: 60, noiseDB: NOISE, keepNoise: false });
  const gC = GRID.xs.map(x => { const c = newCube(GRID.zf.length); E.column(sim, x, GRID.zf, scr0, c); return c; });
  const pC = PROF.xs.map(x => { const c = newCube(PROF.zf.length); E.column(sim, x, PROF.zf, scr0, c); return c; });
  const qC = PROF.qx.map(x => { const c = newCube(PROF.zf.length); E.column(sim, x, PROF.zf, scr0, c); return c; });
  return sets.map(({ k, par }) => {
    const P = Object.assign({}, E.DEF, par), T = {};
    const cols = gC.map(c => E.runColumn(c, P, [k], T, null, null)[k]);
    const pp = pC.map(c => Math.max(E.runColumn(c, P, [k], T, null, [16, 16])[k][16], 0)), pq = qC.map(c => Math.max(E.runColumn(c, P, [k], T, null, [16, 16])[k][16], 0));
    let mx = 0; cols.forEach(col => col.forEach(v => { if (v > mx) mx = v; }));
    const { cy, bg } = roiSplit(cols, GRID.xs, GRID.z0, E.DZ, 0, GRID.nz - 1), rm = E.regionMetrics(cy, bg, mx, [60]), pm = E.profileMetrics(PROF.xs, pp, -3);
    return { k, par, fwhm: pm.fwhm, psl: pm.psl, dip: pairDipM(PROF.qx, pq), dipOld: pairDipOld(PROF.qx, pq), cr: rm.cr, gcnr: rm.gcnr, gcnrD60: rm.disp[0].gcnr, ssnr: rm.ssnr };
  });
}
function displayedAdmire(cls) {
  const sim = E.simulate({ seed: SHOW, ampPoint: 15, ampBright: 60, noiseDB: NOISE, keepNoise: false }), out = [];
  const z0 = GRID.z0, nz = 164, zf = []; for (let i = 0; i <= (nz - 1) * E.SUB; i++) zf.push(z0 + i * E.DZF);
  const cube = newCube(zf.length), pz = PROF.zf.slice(0, 19 * E.SUB + 1), pc = newCube(pz.length);
  for (const cl of cls) {
    const par = Object.assign({}, ADM.RED, { clam: cl }), S = ADM.setup(par), t0 = now(), cols = [];
    for (let c = 0; c < GRID.xs.length; c++) { E.column(sim, GRID.xs[c], zf, scr0, cube); cols.push(Array.from(ADM.line(cube, GRID.xs[c], z0, S, par).subarray(0, GRID.nz))); }
    const prof = { p: [], q: [] };
    for (const [xs, key] of [[PROF.xs, 'p'], [PROF.qx, 'q']]) for (let i = 0; i < xs.length; i++) { E.column(sim, xs[i], pz, scr0, pc); prof[key].push(Math.max(ADM.line(pc, xs[i], PROF.zf[0], S, par, S.hop < S.Lw ? [14, 16] : [14])[16], 0)); }
    let mx = 0; cols.forEach(col => col.forEach(v => { if (v > mx) mx = v; }));
    const { cy, bg } = roiSplit(cols, GRID.xs, GRID.z0, E.DZ, 0, GRID.nz - 1), rm = E.regionMetrics(cy, bg, mx, [60]), pm = E.profileMetrics(PROF.xs, prof.p, -3);
    out.push({ k: 'admire', par: { clam: cl }, fwhm: pm.fwhm, psl: pm.psl, dip: pairDipM(PROF.qx, prof.q), cr: rm.cr, gcnr: rm.gcnr, gcnrD60: rm.disp[0].gcnr, ssnr: rm.ssnr, zc: cy.filter(v => v === 0).length / cy.length, ms: now() - t0 });
    console.log('displayed admire', JSON.stringify(out[out.length - 1]));
  }
  return out;
}

if (part === 'classic' || part === 'all') {
  const ho = heldOutClassic(), lim = FWHM_FACTOR * ho.das.find(r => r.def).fwhm, pick = {};
  Object.keys(ho).forEach(k => { pick[k] = choose(ho[k], lim); const p = pick[k];
    console.log(k.padEnd(7), 'default', JSON.stringify(p.def.par), 'fwhm', p.def.fwhm.toFixed(3), 'gcnr', p.def.gcnr.toFixed(4), '| pick', JSON.stringify(p.pick.par), 'fwhm', p.pick.fwhm.toFixed(3), 'gcnr', p.pick.gcnr.toFixed(4), p.feasible ? '' : 'INFEASIBLE', p.kept ? '(default kept)' : ''); });
  const sets = []; Object.keys(ho).forEach(k => { sets.push({ k, par: pick[k].def.par, role: 'default' }); if (pick[k].pick !== pick[k].def) sets.push({ k, par: pick[k].pick.par, role: 'tuned' }); });
  const t0 = now(), disp = displayedEval(sets); disp.forEach((d, i) => d.role = sets[i].role);
  console.log('displayed eval', ((now() - t0) / 1000).toFixed(1), 's'); disp.forEach(d => console.log(d.k.padEnd(7), d.role.padEnd(7), JSON.stringify(d.par), 'fwhm', d.fwhm.toFixed(3), 'gcnr', d.gcnr.toFixed(4), 'cr', d.cr.toFixed(1), 'ssnr', d.ssnr.toFixed(2), 'dip', d.dip, d.dipOld.toFixed(1)));
  fs.writeFileSync(path.join(OUT, 'classic.json'), JSON.stringify({ HO, SHOW, NOISE, TOL, lim, dasFwhm: ho.das.find(r => r.def).fwhm, ho, pick: Object.fromEntries(Object.entries(pick).map(([k, p]) => [k, { par: p.pick.par, feasible: p.feasible, kept: !!p.kept, same: p.pick === p.def }])), disp }, null, 1));
}
if (part === 'admire' || part === 'all') {
  const ho = heldOutAdmire(), c = fs.existsSync(path.join(OUT, 'classic.json')) ? JSON.parse(fs.readFileSync(path.join(OUT, 'classic.json'), 'utf8')) : { lim: 1.1 * 0.4788 }, p = choose(ho, c.lim);   // lim: 1.1 x DAS FWHM on the held-out data (0.5267 mm, Oct 24 classic run)
  console.log('admire pick', JSON.stringify(p.pick.par), p.feasible, p.kept);
  const cls = [0.0189]; if (p.pick !== p.def) cls.push(p.pick.par.clam);
  const disp = displayedAdmire(cls);
  fs.writeFileSync(path.join(OUT, 'admire.json'), JSON.stringify({ ho, pick: { par: p.pick.par, feasible: p.feasible, kept: !!p.kept, same: p.pick === p.def }, disp }, null, 1));
}
if (part === 'pair') {   // validity test of the pair-dip rule across noise levels, all methods, displayed seed
  for (const nd of [-40, -20, 0, 10, 20, 30]) {
    const sim = E.simulate({ seed: SHOW, ampPoint: 15, ampBright: 60, noiseDB: nd, keepNoise: false });
    const qC = PROF.qx.map(x => { const c = newCube(PROF.zf.length); E.column(sim, x, PROF.zf, scr0, c); return c; }), T = {};
    const r = qC.map(c => E.runColumn(c, E.DEF, GRID.all, T, null, [16, 16]));
    console.log(nd, GRID.all.map(k => { const pq = r.map(o => Math.max(o[k][16], 0)), d = pairDipM(PROF.qx, pq); return k + ' ' + (d === null ? 'n/a' : d.toFixed(1)) + ' (' + pairDipOld(PROF.qx, pq).toFixed(1) + ')'; }).join(' | '));
  }
}
