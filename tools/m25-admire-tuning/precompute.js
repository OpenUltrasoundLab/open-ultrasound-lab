// Offline full-model ADMIRE for the Demo 25.6 grid: noise {-40..30 step 10} dB x aberration {0, 20, 40} ns, shared-recording noise.
// Uses the page's own engine (M16E.simulate/screen/column) and ADMIRE code (ADM) read from m25.html, and the page's own GRID/PROF definitions.
// Usage: node precompute.js <worker k> <of n>   -> appends one JSON line per image line to out/lines_k.jsonl (resumable)
const fs = require('fs'), path = require('path');
const { Lab, loadAll } = require('./load.js'); const { E, ADM } = loadAll();
const page = fs.readFileSync('./m25.html', 'utf8');
const gsrc = page.slice(page.indexOf('  const GRID = { ids:'), page.indexOf('ROI.rows = [8, 61];') + 'ROI.rows = [8, 61];'.length);
const { GRID, PROF } = new Function('E', gsrc + '\nreturn { GRID, PROF, ROI };')(E);
const NOISE = [-40, -30, -20, -10, 0, 10, 20, 30], ABR = [0, 20, 40];
const SET = []; ABR.forEach(ab => NOISE.forEach(nd => SET.push({ nd, ab })));
const par = ADM.FULL, S = ADM.setup(par), N = E.NEL;
const jobs = [];
GRID.xs.forEach((x, i) => jobs.push({ kind: 'grid', i, x }));
PROF.xs.forEach((x, i) => jobs.push({ kind: 'prof', i, x }));
PROF.qx.forEach((x, i) => jobs.push({ kind: 'pair', i, x }));
module.exports = { jobs, SET, GRID, PROF, NOISE, ABR, S };
if (require.main !== module) return;
const k = +process.argv[2], n = +process.argv[3], outDir = path.join(__dirname, 'out'); fs.mkdirSync(outDir, { recursive: true });
const outFile = path.join(outDir, `lines_${k}.jsonl`), done = new Set();
if (fs.existsSync(outFile)) fs.readFileSync(outFile, 'utf8').split('\n').filter(Boolean).forEach(l => { const o = JSON.parse(l); done.add(o.kind + ':' + o.i); });
const sims = NOISE.map(nd => E.simulate({ seed: 7, ampPoint: 15, ampBright: 60, noiseDB: nd, keepNoise: false }));   // same seeds and amplitudes as runAll
const scrs = ABR.map(ab => E.screen(ab, 3, 33));
const tStart = Date.now(); let nDone = 0;
const mine = jobs.filter((j, q) => q % n === k);
for (const job of mine) {
  if (done.has(job.kind + ':' + job.i)) continue;
  const t0 = Date.now();
  let zf, z0, wins, keepRows;
  if (job.kind === 'grid') { z0 = GRID.z0; const nz = 164; zf = []; for (let i = 0; i <= (nz - 1) * E.SUB; i++) zf.push(z0 + i * E.DZF); wins = null; keepRows = [0, 160]; }
  else { z0 = PROF.zf[0]; zf = PROF.zf.slice(0, 19 * E.SUB + 1); wins = S.hop < S.Lw ? [14, 16] : [14]; keepRows = [16, 16]; }
  const cubes = SET.map(s => { const c = { re: new Float32Array(zf.length * N), im: new Float32Array(zf.length * N), nz: zf.length }; E.column(sims[NOISE.indexOf(s.nd)], job.x, zf, scrs[ABR.indexOf(s.ab)], c); return c; });
  const no = Math.floor((zf.length - 1) / E.SUB) + 1, amps = SET.map(() => new Float64Array(no)), das = cubes.map(c => E.bfDAS(c).amp);
  const hop = S.hop || S.Lw, ovl = hop < S.Lw, accs = ovl ? SET.map(() => ({ re: new Float64Array(no), im: new Float64Array(no), n: new Float64Array(no) })) : null;
  const starts = wins || (() => { const a = []; for (let w = ovl ? -(S.Lw - hop) : 0; w < no; w += hop) a.push(w); return a; })();
  const yr = new Float64Array(N), yi = new Float64Array(N); let fits = 0, its = 0, Psum = 0;
  for (const w0 of starts) {
    const sps = cubes.map(c => ADM.spectra(c, w0, S));
    const rec = SET.map(() => ({ R: S.bins.map(() => new Float64Array(N)), I: S.bins.map(() => new Float64Array(N)) }));
    for (let q = 0; q < S.bins.length; q++) {
      const Mo = ADM.model(job.x, z0 + w0 * E.DZ, q, S, par), kb = S.bins[q]; Psum += Mo.P;
      SET.forEach((s, si) => { const sp = sps[si]; for (let e = 0; e < N; e++) { yr[e] = sp.Xr[e * S.Lp + kb]; yi[e] = sp.Xi[e * S.Lp + kb]; }
        const r = ADM.enet(Mo, yr, yi, par); ADM.recon(Mo, r.B, rec[si].R[q], rec[si].I[q]); fits++; its += r.it; });
    }
    SET.forEach((s, si) => ADM.synth(S, rec[si].R, rec[si].I, amps[si], w0, accs ? accs[si] : undefined));
  }
  if (accs) SET.forEach((s, si) => { const A = accs[si]; for (let r = 0; r < no; r++) amps[si][r] = A.n[r] ? Math.hypot(A.re[r], A.im[r]) / (A.n[r] * N) : 0; });
  const o = { kind: job.kind, i: job.i, x: job.x, ms: Date.now() - t0, fits, its, Pmean: Math.round(Psum / (starts.length * S.bins.length)),
    adm: amps.map(a => Array.from(a.slice(keepRows[0], keepRows[1] + 1))), das: das.map(a => Array.from(a.slice(keepRows[0], keepRows[1] + 1))) };
  fs.appendFileSync(outFile, JSON.stringify(o) + '\n'); nDone++;
  console.log(`[w${k}] ${job.kind} ${job.i} x=${job.x} ${(o.ms / 1000).toFixed(1)} s, P~${o.Pmean}, it/fit ${(its / fits).toFixed(1)}; elapsed ${((Date.now() - tStart) / 60000).toFixed(1)} min, ${nDone} lines`);
}
console.log(`[w${k}] finished in ${((Date.now() - tStart) / 60000).toFixed(1)} min`);
