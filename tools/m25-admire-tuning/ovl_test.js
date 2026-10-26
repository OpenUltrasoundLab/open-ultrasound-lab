// STFT overlap test (reduced model, displayed realization): overlap 0 (as released) vs 0.5, at several noise levels
const fs = require('fs');
const { Lab, loadAll } = require('./load.js'); Lab.nextPow2 = n => 1 << Math.ceil(Math.log2(Math.max(2, n)));
const { E, ADM } = loadAll();
const page = fs.readFileSync('./m25.html', 'utf8');
const gsrc = page.slice(page.indexOf('  const GRID = { ids:'), page.indexOf('ROI.rows = [8, 61];') + 'ROI.rows = [8, 61];'.length);
const { GRID, PROF } = new Function('E', gsrc + '\nreturn { GRID, PROF, ROI };')(E);
const N = E.NEL, newCube = nz => ({ re: new Float32Array(nz * N), im: new Float32Array(nz * N), nz }), scr0 = E.screen(0, 3, 33);
function roiSplit(cols) { const rc = E.scene.roiCyst, rb = E.scene.roiBg, cy = [], bg = [];
  for (let c = 0; c < GRID.xs.length; c++) { const x = GRID.xs[c]; for (let r = 0; r < GRID.nz; r++) { const z = Math.round((GRID.z0 + r * E.DZ) * 1e4) / 1e4, v = Math.max(cols[c][r], 0);
    if ((x - rc[0]) ** 2 + (z - rc[1]) ** 2 < rc[2] ** 2) cy.push(v); if (x >= rb[0] - 1e-9 && x <= rb[1] + 1e-9 && z >= rb[2] - 1e-9 && z <= rb[3] + 1e-9) bg.push(v); } } return { cy, bg }; }
const noises = (process.argv[2] || '-20,0').split(',').map(Number), clam = 0.0189 / 4;
for (const nd of noises) {
  const sim = E.simulate({ seed: 7, ampPoint: 15, ampBright: 60, noiseDB: nd, keepNoise: false });
  const z0 = GRID.z0, nz = 164, zf = []; for (let i = 0; i <= (nz - 1) * E.SUB; i++) zf.push(z0 + i * E.DZF);
  const cube = newCube(zf.length), pz = PROF.zf.slice(0, 19 * E.SUB + 1), pc = newCube(pz.length);
  const dcols = []; for (const x of GRID.xs) { E.column(sim, x, GRID.zf, scr0, cube); }
  for (const ov of [0, 0.5]) {
    const par = Object.assign({}, ADM.RED, { clam, overlap: ov }), S = ADM.setup(par), t0 = Date.now(), cols = [];
    for (const x of GRID.xs) { E.column(sim, x, zf, scr0, cube); cols.push(Array.from(ADM.line(cube, x, z0, S, par).subarray(0, GRID.nz))); }
    const pw = ov ? [14, 16] : [14], prof = [];
    for (const x of PROF.xs) { E.column(sim, x, pz, scr0, pc); prof.push(Math.max(ADM.line(pc, x, PROF.zf[0], S, par, pw)[16], 0)); }
    let mx = 0; cols.forEach(col => col.forEach(v => { if (v > mx) mx = v; }));
    const { cy, bg } = roiSplit(cols), rm = E.regionMetrics(cy, bg, mx, [60]), pm = E.profileMetrics(PROF.xs, prof, -3);
    const zbg = bg.filter(v => v === 0).length / bg.length;
    console.log(JSON.stringify({ nd, overlap: ov, gcnr: +rm.gcnr.toFixed(3), cr: +rm.cr.toFixed(1), ssnr: +rm.ssnr.toFixed(2), fwhm: +pm.fwhm.toFixed(3), zeroCyst: +(cy.filter(v => v === 0).length / cy.length).toFixed(2), zeroBg: +zbg.toFixed(3), s: (Date.now() - t0) / 1000 }));
  }
}
