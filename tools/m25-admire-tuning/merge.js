// Merge out/lines_*.jsonl -> admire_grid.json: per setting an 8-bit dB image (81 x 161, same display mapping as the page: 20 log10(v / max),
// stored as q = round(255 (dB + 100) / 100), q = 0 for <= -100 dB or exact zero) and metrics from the page's own metric functions
const fs = require('fs'), path = require('path');
const { Lab, loadAll } = require('./load.js'); const { E, ADM } = loadAll();
const P = require('./precompute.js');
const page = fs.readFileSync('./m25.html', 'utf8');
const grabFn = name => { const i = page.indexOf('  function ' + name + '('); let j = page.indexOf('\n  function ', i + 10); const k = page.indexOf('\n  const ', i + 10); if (k > 0 && k < j) j = k; return page.slice(i, j); };
const roiSplit = new Function('E', grabFn('roiSplit') + '\nreturn roiSplit;')(E);
const pairDip = new Function(grabFn('pairDip') + '\nreturn pairDip;')();
const { GRID, PROF, SET } = P, NX = GRID.xs.length, NZ = GRID.nz, DRS = []; for (let d = 20; d <= 100; d += 5) DRS.push(d);
const lines = {}; let nl = 0, msTot = 0, fitsTot = 0, msGrid = 0, msProf = 0;
fs.readdirSync(path.join(__dirname, 'out')).filter(f => /^lines_\d+\.jsonl$/.test(f)).forEach(f => fs.readFileSync(path.join(__dirname, 'out', f), 'utf8').split('\n').filter(Boolean).forEach(l => {
  const o = JSON.parse(l); lines[o.kind + ':' + o.i] = o; nl++; msTot += o.ms; fitsTot += o.fits; if (o.kind === 'grid') msGrid += o.ms; else msProf += o.ms; }));
const have = k => !!lines[k];
const nGrid = GRID.xs.filter((x, i) => have('grid:' + i)).length, nProf = PROF.xs.filter((x, i) => have('prof:' + i)).length, nPair = PROF.qx.filter((x, i) => have('pair:' + i)).length;
console.log(`lines: grid ${nGrid}/${NX}, prof ${nProf}/${PROF.xs.length}, pair ${nPair}/${PROF.qx.length}; CPU ${(msTot / 3.6e6).toFixed(2)} h (grid ${(msGrid / 3.6e6).toFixed(2)} h, profiles ${(msProf / 3.6e6).toFixed(2)} h), ${fitsTot} fits`);
function metrics(cols, prof, q) {
  let mx = 0; cols.forEach(c => c.forEach(v => { if (v > mx) mx = v; }));
  const { cy, bg } = roiSplit(cols, GRID.xs, GRID.z0, E.DZ, 0, NZ - 1), rm = E.regionMetrics(cy, bg, mx, DRS);
  const m = { cr: rm.cr, gcnr: rm.gcnr, ssnr: rm.ssnr, zc: cy.filter(v => v === 0).length / cy.length, zb: bg.filter(v => v === 0).length / bg.length,
    cnrD: rm.disp.map(d => d.cnr), crD: rm.disp.map(d => d.cr), gcnrD: rm.disp.map(d => d.gcnr) };
  if (prof && q) { const pm = E.profileMetrics(PROF.xs, prof, -3); m.fwhm = pm.fwhm; m.psl = pm.psl; m.dip = pairDip(PROF.qx, q); }
  return m;
}
const out = { v: 1, noise: P.NOISE, ab: P.ABR, nx: NX, nz: NZ, DRS, complete: nGrid === NX && nProf === PROF.xs.length && nPair === PROF.qx.length,
  cpuHours: msTot / 3.6e6, fits: fitsTot, settings: [] };
const imgs = [];
SET.forEach((s, si) => {
  const colsA = [], colsD = [];
  for (let i = 0; i < NX; i++) { const o = lines['grid:' + i]; colsA.push(o ? o.adm[si] : new Array(NZ).fill(0)); colsD.push(o ? o.das[si] : new Array(NZ).fill(0)); }
  const pa = PROF.xs.map((x, i) => lines['prof:' + i] ? Math.max(lines['prof:' + i].adm[si][0], 0) : 0), pq = PROF.qx.map((x, i) => lines['pair:' + i] ? Math.max(lines['pair:' + i].adm[si][0], 0) : 0);
  const da = PROF.xs.map((x, i) => lines['prof:' + i] ? lines['prof:' + i].das[si][0] : 0), dq = PROF.qx.map((x, i) => lines['pair:' + i] ? lines['pair:' + i].das[si][0] : 0);
  const mA = metrics(colsA, nProf === PROF.xs.length ? pa : null, nPair === PROF.qx.length ? pq : null), mD = metrics(colsD, nProf === PROF.xs.length ? da : null, nPair === PROF.qx.length ? dq : null);
  let mx = 0; colsA.forEach(c => c.forEach(v => { if (v > mx) mx = v; }));
  const img = new Array(NX * NZ); for (let r = 0; r < NZ; r++) for (let i = 0; i < NX; i++) { const v = colsA[i][r]; const db = v > 0 && mx > 0 ? 20 * Math.log10(v / mx) : -1e9; img[r * NX + i] = Math.max(0, Math.min(255, Math.round(255 * (db + 100) / 100))); }
  imgs.push(img);
  const r4 = v => v === undefined ? null : Math.round(v * 1e4) / 1e4;
  const pack = m => ({ cr: r4(m.cr), gcnr: r4(m.gcnr), ssnr: r4(m.ssnr), fwhm: r4(m.fwhm), psl: r4(m.psl), dip: r4(m.dip), zc: r4(m.zc), zb: r4(m.zb), cnrD: m.cnrD.map(r4), gcnrD: m.gcnrD.map(r4), crD: m.crD.map(r4) });
  out.settings.push({ nd: s.nd, ab: s.ab, adm: pack(mA), das: pack(mD) });
});
// mean per-image times (one setting = 1/24 of the shared-model work is not meaningful; report model-shared totals)
out.msPerLineGrid = msGrid / Math.max(1, nGrid); out.msProfAll = msProf;
fs.writeFileSync(path.join(__dirname, 'admire_grid.json'), JSON.stringify(out));
fs.writeFileSync(path.join(__dirname, 'admire_imgs.json'), JSON.stringify(imgs));
const t = s => `${s.nd}dB/${s.ab}ns`;
out.settings.forEach(s => console.log(t(s).padEnd(10), 'ADMIRE CR', s.adm.cr, 'gCNR', s.adm.gcnr, 'sSNR', s.adm.ssnr, 'FWHM', s.adm.fwhm, 'dip', s.adm.dip, 'zeroCyst', s.adm.zc, ' | DAS CR', s.das.cr, 'gCNR', s.das.gcnr, 'sSNR', s.das.ssnr, 'FWHM', s.das.fwhm));
