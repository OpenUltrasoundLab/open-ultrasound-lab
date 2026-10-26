// Validation (a): on-axis point target + strong off-axis point clutter (+ noise), line x = 0; full and reduced models (page code)
const { Lab, loadAll } = require('./load.js'); const { E, ADM } = loadAll();
const TAU = 2 * Math.PI, NT = E.NT;
function simPts(pts, noiseRel, seed) {   // pts: [x, z, amp]; noise: complex Gaussian filtered by the pulse, power = noiseRel x (unit point echo peak power)
  const A_ = 1 / (2 * E.SIGT ** 2), K = Math.ceil(4 * E.SIGT * E.FS), rng = Lab.rng(seed || 5);
  const taps = []; let e2 = 0; for (let k = -K; k <= K; k++) { const g = Math.exp(-A_ * (k / E.FS) ** 2); taps.push(g); e2 += g * g; }
  const sig = Math.sqrt(noiseRel / e2);
  return { data: E.ANG.map(th => { const re = new Float32Array(64 * NT), im = new Float32Array(64 * NT);
    for (const [xs, zs, a] of pts) { const txd = (xs * Math.sin(th) + zs * Math.cos(th)) / E.C;
      for (let e = 0; e < 64; e++) { const dx = xs - E.XE[e], tau = txd + Math.sqrt(dx * dx + zs * zs) / E.C, ph = -TAU * E.F0 * tau, n0 = (tau - E.T0) * E.FS;
        for (let n = Math.max(0, Math.ceil(n0 - K)); n <= Math.min(NT - 1, Math.floor(n0 + K)); n++) { const u = (n - n0) / E.FS, g = a * Math.exp(-A_ * u * u); re[e * NT + n] += g * Math.cos(ph); im[e * NT + n] += g * Math.sin(ph); } } }
    if (noiseRel > 0) { const wr = new Float64Array(NT + 2 * K), wi = new Float64Array(NT + 2 * K);
      for (let e = 0; e < 64; e++) { for (let n = 0; n < wr.length; n++) { const z = rng.cg(); wr[n] = z.re * sig; wi[n] = z.im * sig; }
        for (let n = 0; n < NT; n++) { let sr = 0, si = 0; for (let k = 0; k <= 2 * K; k++) { sr += taps[k] * wr[n + k]; si += taps[k] * wi[n + k]; } re[e * NT + n] += sr; im[e * NT + n] += si; } } }
    return { re, im }; }) };
}
const target = [[0, 30, 1]], clutter = [[1.5, 30.02, 3], [2.5, 29.95, 3], [-2.0, 30.05, 3], [-3.5, 29.97, 3], [4.5, 30.1, 3]];
const scenes = { C: simPts(target, 0), B: simPts(clutter, 1e-4, 9), A: simPts(target.concat(clutter), 1e-4, 9) };
const scr = E.screen(0, 3, 33), z0 = 28.6, zf = []; for (let i = 0; i <= (40 - 1) * E.SUB; i++) zf.push(z0 + i * E.DZF);   // windows from 28.6 mm: 29.8-30.025 contains z = 30
const rowsZ = Array.from({ length: 40 }, (_, r) => z0 + r * E.DZ), it = rowsZ.findIndex(z => Math.abs(z - 30) < 0.04);
const res = {};
for (const [name, par] of [['full', ADM.FULL], ['reduced', ADM.RED]]) { const S = ADM.setup(par); res[name] = {};
  for (const [k, sim] of Object.entries(scenes)) { const cube = { re: new Float32Array(zf.length * 64), im: new Float32Array(zf.length * 64), nz: zf.length }; E.column(sim, 0, zf, scr, cube);
    const das = E.bfDAS(cube).amp, adm = ADM.line(cube, 0, z0, S, par); res[name][k] = { das: Array.from(das), adm: Array.from(adm) }; } }
const db = v => (20 * Math.log10(v)).toFixed(1);
for (const name of ['full', 'reduced']) { const R = res[name], dT = R.C.das[it], aT = R.C.adm[it];
  const away = rowsZ.map((z, r) => Math.abs(z - 30) > 0.45 ? r : -1).filter(r => r >= 0);
  const mean2 = (a, rows) => Math.sqrt(rows.reduce((s, r) => s + a[r] * a[r], 0) / rows.length);
  console.log(`[${name}] target only: ADMIRE/DAS at the target = ${db(aT / dT)} dB`);
  console.log(`[${name}] clutter only at the target pixel: DAS ${db(R.B.das[it] / dT)} dB, ADMIRE ${R.B.adm[it] > 0 ? db(R.B.adm[it] / aT) + ' dB' : 'exactly 0'} (re target-only peak of the same method)`);
  console.log(`[${name}] clutter only, rms over rows > 0.45 mm from 30 mm: DAS ${db(mean2(R.B.das, away) / dT)} dB, ADMIRE ${mean2(R.B.adm, away) > 0 ? db(mean2(R.B.adm, away) / aT) + ' dB' : 'exactly 0'}`);
  console.log(`[${name}] target + clutter: DAS ${db(R.A.das[it] / dT)} dB, ADMIRE ${db(R.A.adm[it] / aT)} dB re own target-only value; ADMIRE zero rows elsewhere: ${away.filter(r => R.A.adm[r] === 0).length}/${away.length}`);
}
// lateral profile through the target alone (full model): FWHM of DAS and ADMIRE
const xs = []; for (let x = -1.2; x <= 1.2001; x += 0.02) xs.push(+x.toFixed(3));
const zw = []; for (let i = 0; i <= 3 * E.SUB; i++) zw.push(29.8 + i * E.DZF);
for (const [name, par] of [['full', ADM.FULL], ['reduced', ADM.RED]]) { const S = ADM.setup(par), pd = [], pa = [];
  for (const x of xs) { const cube = { re: new Float32Array(zw.length * 64), im: new Float32Array(zw.length * 64), nz: zw.length }; E.column(scenes.C, x, zw, scr, cube);
    pd.push(E.bfDAS(cube).amp[3]); pa.push(ADM.line(cube, x, 29.8, S, par)[3]); }   // row 3 = 30.025 mm
  const pmD = E.profileMetrics(xs, pd, 0), pmA = E.profileMetrics(xs, pa, 0);
  console.log(`[${name}] lateral profile at z = 30.025 mm, target alone: FWHM DAS ${pmD.fwhm.toFixed(3)} mm, ADMIRE ${pmA.fwhm.toFixed(3)} mm; ADMIRE zero beyond |x| = ${(() => { let m = 0; xs.forEach((x, i) => { if (pa[i] > 0) m = Math.max(m, Math.abs(x)); }); return m.toFixed(2); })()} mm`);
}
