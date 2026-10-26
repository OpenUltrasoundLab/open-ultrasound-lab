// Per-bin wavenumber calibration: maximize mean (over isolated point scatterers) of the best model-data correlation over d
const { Lab, loadAll } = require('./load.js'); const { E, ADM } = loadAll();
const S = ADM.setup(ADM.FULL), TAU = 2 * Math.PI, NT = E.NT;
function simPoint(xs, zs) { const A_ = 1 / (2 * E.SIGT ** 2), K = Math.ceil(4 * E.SIGT * E.FS);
  return { data: E.ANG.map(th => { const re = new Float32Array(64 * NT), im = new Float32Array(64 * NT), txd = (xs * Math.sin(th) + zs * Math.cos(th)) / E.C;
    for (let e = 0; e < 64; e++) { const dx = xs - E.XE[e], tau = txd + Math.sqrt(dx * dx + zs * zs) / E.C, ph = -TAU * E.F0 * tau, n0 = (tau - E.T0) * E.FS;
      for (let n = Math.max(0, Math.ceil(n0 - K)); n <= Math.min(NT - 1, Math.floor(n0 + K)); n++) { const u = (n - n0) / E.FS, g = Math.exp(-A_ * u * u); re[e * NT + n] += g * Math.cos(ph); im[e * NT + n] += g * Math.sin(ph); } }
    return { re, im }; }) }; }
const scr = E.screen(0, 3, 33), cases = [];
for (const [xcol, zfirst] of [[0, 29.8875], [-4.5, 23.8875]]) for (const dx of [0, 0.25, 0.5, 1, 2, 3, 5]) for (const dz of [-0.1, 0, 0.1]) {
  const zc = zfirst + 1.5 * E.DZ, xs = xcol + dx, zs = zc + dz, sim = simPoint(xs, zs), zf = []; for (let i = 0; i < 4 * E.SUB; i++) zf.push(zfirst + i * E.DZF);
  const cube = { re: new Float32Array(zf.length * 64), im: new Float32Array(zf.length * 64), nz: zf.length }; E.column(sim, xcol, zf, scr, cube);
  const sp = ADM.spectra(cube, 0, S); cases.push({ xcol, zfirst, zc, xs, zs, sp }); }
const erf = ADM.erf, h = S.Lw / ADM.FSA * E.C / 2, st = E.C * E.SIGT;
function corr(cs, q, cal) { const k = S.bins[q], kw = TAU * S.freqs[q] * cal / E.C, yr = [], yi = []; let ny = 0;
  for (let e = 0; e < 64; e++) { yr.push(cs.sp.Xr[e * S.Lp + k]); yi.push(cs.sp.Xi[e * S.Lp + k]); ny += yr[e] ** 2 + yi[e] ** 2; } ny = Math.sqrt(ny);
  const z2a = 2 * cs.zfirst, z2b = 2 * (cs.zfirst + 3 * E.DZ); let best = 0;
  for (let d = -2; d <= 2; d += 0.004) { let cr = 0, ci = 0, nt = 0, mx = 0; const a = [], zd = [];
    for (let e = 0; e < 64; e++) { const u = E.XE[e] - cs.xcol, R = Math.hypot(u - (cs.xs - cs.xcol), cs.zs), d0 = R + d + 2 * cs.zc - cs.zs, z = (d0 * d0 - u * u) / d0; zd.push(z);
      let g = 0; if (!(z - h > z2b || z + h < z2a)) { const p0 = Math.max(z - h, z2a), p1 = Math.min(z + h, z2b); g = erf((z - p0) / st) - erf((z - p1) / st); } a.push(g); if (g > mx) mx = g; }
    if (!(mx > 0)) continue;
    for (let e = 0; e < 64; e++) { const A = Math.sqrt(a[e] / mx), ph = -kw * zd[e], tr = A * Math.cos(ph), ti = A * Math.sin(ph); nt += A * A; cr += tr * yr[e] + ti * yi[e]; ci += tr * yi[e] - ti * yr[e]; }
    const c = Math.hypot(cr, ci) / ny / Math.sqrt(nt); if (c > best) best = c; }
  return best; }
const res = [];
for (let q = 0; q < S.bins.length; q++) { let bc = 0, bm = -1; const tab = [];
  for (let cal = 0.80; cal <= 1.2001; cal += 0.01) { const m = cases.reduce((s, cs) => s + corr(cs, q, cal), 0) / cases.length; tab.push(cal.toFixed(2) + ':' + m.toFixed(4)); if (m > bm) { bm = m; bc = cal; } }
  const m1 = cases.reduce((s, cs) => s + corr(cs, q, 1), 0) / cases.length;
  const inROI = cases.filter(c => Math.abs(c.xs - c.xcol) <= 0.25), off = cases.filter(c => Math.abs(c.xs - c.xcol) >= 1);
  res.push({ f: S.freqs[q], cal: +bc.toFixed(2), meanCorrCal: bm, meanCorr1: m1, inROI1: inROI.reduce((s, c) => s + corr(c, q, 1), 0) / inROI.length, inROIcal: inROI.reduce((s, c) => s + corr(c, q, bc), 0) / inROI.length, off1: off.reduce((s, c) => s + corr(c, q, 1), 0) / off.length, offcal: off.reduce((s, c) => s + corr(c, q, bc), 0) / off.length });
  console.log(tab.join(' ')); }
console.log(JSON.stringify(res, null, 1));
