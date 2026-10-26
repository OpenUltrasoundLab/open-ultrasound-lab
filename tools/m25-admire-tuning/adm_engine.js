/* ===== ADMIRE (aperture domain model image reconstruction; Byram et al. 2015), port of the authors' released code (VU-BEAM-Lab/ADMIRE),
   adapted to this demo's delayed, angle-compounded analytic aperture data (no DOM; shared by the page and the offline precompute) ===== */
const ADM = (function (E) {
  const C = E.C, TAU = 2 * Math.PI, NEL = E.NEL, XE = E.XE, PITCH = E.PITCH, DZ = E.DZ, SUB = E.SUB;
  const FSA = C / (2 * DZ);   // equivalent two-way sampling of the image axial grid: 10.27 MHz (complex samples)
  const ERF_T = (() => { const n = 60001, t = new Float64Array(n); for (let i = 0; i < n; i++) { const x = i * 1e-4, q = 1 / (1 + 0.3275911 * x);   // A&S 7.1.26 (|err| < 1.5e-7), tabulated on [0, 6]
    t[i] = 1 - (((((1.061405429 * q - 1.453152027) * q) + 1.421413741) * q - 0.284496736) * q + 0.254829592) * q * Math.exp(-x * x); } return t; })();
  const erf = x => { const a = x < 0 ? -x : x; if (a >= 6) return x < 0 ? -1 : 1; const p = a * 1e4, i = p | 0, v = ERF_T[i] + (p - i) * (ERF_T[i + 1] - ERF_T[i]); return x < 0 ? -v : v; };
  const colon = (a, s, b) => { const n = Math.floor((b - a) / s + 1e-9) + 1, o = []; for (let i = 0; i < n; i++) o.push(a + i * s); return o; };   // MATLAB a:s:b
  /* FULL = released defaults (lengths in mm): ROI and outer model spaces, alpha, c_lambda; outer depths = README's off-axis-only option.
     kcal = per-bin wavenumber calibration (the release uses empirical per-probe factors; ours are fitted to this simulator, see the notes). */
  const FULL = { name: 'full', pulseMult: 1, span: 1.2, bins: null, roiX: [-0.5, 0.0179, 0.5], roiZ: [-0.5, 0.1430, 0.5], roiD: [-0.8, 0.0485, 0.4],
    outX: [-1, 1.4228, 1], outZ: [1, -0.5, 0.7114, 1, 0.5, 0], outD: [-8, 0.1211, 3.2], latOff: 1, ell1: 0, ell2: 0.5, alpha: 0.9, clam: 0.0189, tol: 1e-5, maxit: 100000,
    kcal: [1.10, 0.95, 0.88] };
  function setup(par) {
    const fwhm = Math.sqrt(8 * Math.LN2) * E.SIGT, Lw = Math.ceil(par.pulseMult * fwhm * FSA), Lp = 2 * Lw;
    const lo = E.F0 * (1 - 0.5 * E.BW * par.span), hi = E.F0 * (1 + 0.5 * E.BW * par.span); let bins = [];
    for (let k = 0; k < Lp; k++) { const f = k * FSA / Lp; if (f >= lo - 1e-9 && f <= hi + 1e-9) bins.push(k); }
    const all = bins.slice(); if (par.bins) bins = par.bins.map(i => all[i]);
    return { Lw, Lp, bins, all, freqs: bins.map(k => k * FSA / Lp), kcal: bins.map(k => par.kcal[all.indexOf(k)]), fwhm, FSA };
  }
  /* aperture signal of a scatterer at (xp, zp) (xp relative to the image line) with extra path d, for a window [z2a, z2b] (two-way distance);
     u: element positions relative to the line; Rx: |element - scatterer| (precomputed); writes amplitude A_FT and zd; returns max power or 0 */
  function sig(u, Rx, zp, d, zc, z2a, z2b, h, st, ga, zd) {
    let mx = 0;
    for (let e = 0; e < NEL; e++) {
      const d0 = Rx[e] + d + 2 * zc - zp, z = (d0 * d0 - u[e] * u[e]) / d0; zd[e] = z;   // z = 2 z', z' = perceived depth (dynamic receive focus)
      if (z - h > z2b || z + h < z2a) { ga[e] = 0; continue; }
      const pmin = Math.max(z - h, z2a), pmax = Math.min(z + h, z2b), g = (erf((z - pmin) / st) - erf((z - pmax) / st));
      ga[e] = g; if (g > mx) mx = g;
    }
    return mx;
  }
  /* model matrix for one line (xcol), one STFT window (first row depth zfirst), one bin index bi: complex unit-norm columns, ROI first */
  function model(xcol, zfirst, bi, S, par) {
    const Lw = S.Lw, zc = zfirst + (Lw - 1) / 2 * DZ, z2a = 2 * zfirst, z2b = 2 * (zfirst + (Lw - 1) * DZ), h = Lw / FSA * C / 2;
    const f = S.freqs[bi], kw = TAU * f * S.kcal[bi] / C, lam = TAU / kw, st = C * E.SIGT;
    const u = new Float64Array(NEL); for (let e = 0; e < NEL; e++) u[e] = XE[e] - xcol;
    const ga = new Float64Array(NEL), zd = new Float64Array(NEL), Rx = new Float64Array(NEL);
    // lateral resolution scale from the FFT (4096 points) of |on-axis model|, as in the release
    for (let e = 0; e < NEL; e++) Rx[e] = Math.sqrt(u[e] * u[e] + zc * zc);
    const M = 1 << (2 * Math.floor(Math.log2(NEL))), fr = new Float64Array(M), fi = new Float64Array(M), m0 = sig(u, Rx, zc, 0, zc, z2a, z2b, h, st, ga, zd);
    for (let e = 0; e < NEL; e++) fr[e] = Math.sqrt(ga[e] / m0); Lab.fft(fr, fi, false);
    let mx = 0; for (let k = 0; k < M; k++) { fr[k] = Math.hypot(fr[k], fi[k]); if (fr[k] > mx) mx = fr[k]; }
    let fw = 0; for (let k = 0; k < M; k++) if (fr[k] > 0.5 * mx) fw++;
    const resLat = fw / PITCH / M * lam * zc, resAx = 2 * resLat;
    let cap = 4096, P = 0, xr = new Float32Array(cap * NEL), xi = new Float32Array(cap * NEL);
    const space = (xs, zs, dlo, dstep, dhi, keep) => {
      for (const xx of xs) for (const zz of zs) {
        if (zz <= 0 || !keep(xx, zz)) continue;
        for (let e = 0; e < NEL; e++) Rx[e] = Math.sqrt((u[e] - xx) * (u[e] - xx) + zz * zz);
        // d-range where the pulse overlaps the window at some element (zd increases with d): skip the all-zero predictors cheaply
        let a = Infinity, b = -Infinity;
        for (let e = 0; e < NEL; e++) { const base = Rx[e] + 2 * zc - zz, v1 = z2a - h, v2 = z2b + h, uu = 4 * u[e] * u[e];
          const lo = (v1 + Math.sqrt(v1 * v1 + uu)) / 2 - base, hi = (v2 + Math.sqrt(v2 * v2 + uu)) / 2 - base; if (lo < a) a = lo; if (hi > b) b = hi; }
        const n = Math.floor((dhi - dlo) / dstep + 1e-9) + 1, i0 = Math.max(0, Math.ceil((a - dlo) / dstep - 1e-9)), i1 = Math.min(n - 1, Math.floor((b - dlo) / dstep + 1e-9));
        for (let i = i0; i <= i1; i++) {
          const dd = dlo + i * dstep, m = sig(u, Rx, zz, dd, zc, z2a, z2b, h, st, ga, zd); if (!(m > 0)) continue;
          if (P === cap) { cap *= 2; const r2 = new Float32Array(cap * NEL), i2 = new Float32Array(cap * NEL); r2.set(xr); i2.set(xi); xr = r2; xi = i2; }
          const o = P * NEL; let s = 0;
          for (let e = 0; e < NEL; e++) { const A = Math.sqrt(ga[e] / m); if (!A) { xr[o + e] = 0; xi[o + e] = 0; continue; } const ph = -kw * zd[e]; xr[o + e] = A * Math.cos(ph); xi[o + e] = A * Math.sin(ph); s += A * A; }
          s = 1 / Math.sqrt(s); for (let e = 0; e < NEL; e++) { xr[o + e] *= s; xi[o + e] *= s; }
          P++;
        }
      }
    };
    const ex = par.ell1, ox = par.ell2;
    space(colon(par.roiX[0] * resLat, par.roiX[1] * resLat, par.roiX[2] * resLat), colon(zc + par.roiZ[0] * resAx, par.roiZ[1] * resAx, zc + par.roiZ[2] * resAx),
      par.roiD[0], par.roiD[1] * lam, par.roiD[2], (x, z) => (x / (resLat + ex)) ** 2 + ((z - zc) / (resAx + ex)) ** 2 <= 1);
    const nROI = P;
    const L = PITCH * NEL / 2 + par.latOff, oz = par.outZ;
    space(colon(par.outX[0] * L, par.outX[1] * resLat, par.outX[2] * L), colon(oz[0] * zc + oz[1] * resAx, oz[2] * resAx, oz[3] * zc + oz[4] * resAx + oz[5] * zc),
      par.outD[0], par.outD[1] * lam, par.outD[2], (x, z) => (x / (resLat + ox)) ** 2 + ((z - zc) / (resAx + ox)) ** 2 > 1);
    xr = xr.subarray(0, P * NEL); xi = xi.subarray(0, P * NEL);
    return { xr, xi, P, nROI, resLat, resAx, zc, f };
  }
  /* elastic net by cyclic coordinate descent on the real-ified complex model (as ccd_double_precision.c):
     min (1/2No)||y - X b||^2 + lam (alpha ||b||_1 + (1 - alpha)||b||^2 / 2), No = 2 NEL, unit-norm columns, y standardized, lam = c_lam rms(y);
     sweeps over the active set between full sweeps (same fixed point; stop when a full sweep changes no coefficient by more than tol);
     B: optional warm start (unstandardized), length 2P: [Re, Im] coefficient of each complex predictor */
  function enet(Mo, yr, yi, par, B) {
    const P = Mo.P, xr = Mo.xr, xi = Mo.xi, No = 2 * NEL; let ss = 0, m = 0;
    for (let e = 0; e < NEL; e++) { ss += yr[e] * yr[e] + yi[e] * yi[e]; m += yr[e] + yi[e]; }
    m /= No; const rms = Math.sqrt(ss / No), sd = Math.sqrt(Math.max(ss / No - m * m, 0));
    if (!B || B.length !== 2 * P) B = new Float64Array(2 * P);
    if (!(sd > 0)) { B.fill(0); return { B, it: 0 }; }
    const lam = par.clam * rms / sd, g = lam * par.alpha, den = 1 / No + lam * (1 - par.alpha), tol = par.tol;
    const rr = new Float64Array(NEL), ri = new Float64Array(NEL);
    for (let e = 0; e < NEL; e++) { rr[e] = yr[e] / sd; ri[e] = yi[e] / sd; }
    for (let j = 0; j < P; j++) { const a = (B[2 * j] /= sd), b = (B[2 * j + 1] /= sd); if (a || b) { const o = j * NEL; for (let e = 0; e < NEL; e++) { rr[e] -= xr[o + e] * a - xi[o + e] * b; ri[e] -= xi[o + e] * a + xr[o + e] * b; } } }
    const upd = j => {
      const o = j * NEL; let mc = 0, a = B[2 * j], p = 0;
      for (let e = 0; e < NEL; e++) p += xr[o + e] * rr[e] + xi[o + e] * ri[e];
      p = (p + a) / No; let na = p > g ? (p - g) / den : p < -g ? (p + g) / den : 0;
      if (na !== a) { const d = na - a; for (let e = 0; e < NEL; e++) { rr[e] -= xr[o + e] * d; ri[e] -= xi[o + e] * d; } B[2 * j] = na; mc = d * d / No; }
      let b = B[2 * j + 1]; p = 0;
      for (let e = 0; e < NEL; e++) p += xr[o + e] * ri[e] - xi[o + e] * rr[e];
      p = (p + b) / No; let nb = p > g ? (p - g) / den : p < -g ? (p + g) / den : 0;
      if (nb !== b) { const d = nb - b; for (let e = 0; e < NEL; e++) { rr[e] += xi[o + e] * d; ri[e] -= xr[o + e] * d; } B[2 * j + 1] = nb; const c = d * d / No; if (c > mc) mc = c; }
      return mc;
    };
    let it = 0;
    while (it < par.maxit) {
      let mc = 0; for (let j = 0; j < P; j++) { const c = upd(j); if (c > mc) mc = c; } it++;
      if (mc < tol) break;
      const act = []; for (let j = 0; j < P; j++) if (B[2 * j] || B[2 * j + 1]) act.push(j);
      while (it < par.maxit) { let m2 = 0; for (let q = 0; q < act.length; q++) { const c = upd(act[q]); if (c > m2) m2 = c; } it++; if (m2 < tol) break; }
    }
    for (let j = 0; j < 2 * P; j++) B[j] *= sd;
    return { B, it };
  }
  /* decluttered aperture vector at one bin: X_ROI b_ROI; also the outer (clutter) part if wanted */
  function recon(Mo, B, outR, outI, clR, clI) {
    outR.fill(0); outI.fill(0); if (clR) { clR.fill(0); clI.fill(0); }
    for (let j = 0; j < Mo.P; j++) { const a = B[2 * j], b = B[2 * j + 1]; if (!a && !b) continue; const o = j * NEL, roi = j < Mo.nROI;
      if (!roi && !clR) continue; const R = roi ? outR : clR, I = roi ? outI : clI;
      for (let e = 0; e < NEL; e++) { R[e] += Mo.xr[o + e] * a - Mo.xi[o + e] * b; I[e] += Mo.xi[o + e] * a + Mo.xr[o + e] * b; } }
  }
  /* STFT of one window (coarse rows w0..w0+Lw-1 of a cube): returns per element the Lp-point spectrum, layout [e][k] */
  function spectra(cube, w0, S) {
    const Lp = S.Lp, Xr = new Float64Array(NEL * Lp), Xi = new Float64Array(NEL * Lp), tr = new Float64Array(Lp), ti = new Float64Array(Lp), no = Math.floor((cube.nz - 1) / SUB) + 1;
    for (let e = 0; e < NEL; e++) { tr.fill(0); ti.fill(0);
      for (let n = 0; n < S.Lw; n++) { const r = w0 + n; if (r >= no) break; const b = r * SUB * NEL + e; tr[n] = cube.re[b]; ti[n] = cube.im[b]; }
      Lab.fft(tr, ti, false); Xr.set(tr, e * Lp); Xi.set(ti, e * Lp); }
    return { Xr, Xi };
  }
  /* inverse STFT of the reconstructed bins (all other bins zero), first Lw samples, then DAS over the aperture: |sum_e| / NEL per row */
  function synth(S, recR, recI, ampOut, w0) {   // recR/recI: per selected bin, Float64Array(NEL)
    const Lp = S.Lp, tr = new Float64Array(Lp), ti = new Float64Array(Lp), sr = new Float64Array(S.Lw), si = new Float64Array(S.Lw);
    for (let e = 0; e < NEL; e++) { tr.fill(0); ti.fill(0); S.bins.forEach((k, q) => { tr[k] = recR[q][e]; ti[k] = recI[q][e]; }); Lab.fft(tr, ti, true);
      for (let n = 0; n < S.Lw; n++) { sr[n] += tr[n]; si[n] += ti[n]; } }
    for (let n = 0; n < S.Lw; n++) if (w0 + n < ampOut.length) ampOut[w0 + n] = Math.hypot(sr[n], si[n]) / NEL;
  }
  /* ADMIRE on one image line: cube = delayed aperture data at the line (fine rows; coarse rows r*SUB are used), zf0 = depth of row 0.
     Windows tile the line from row 0 (no overlap, as released); win: optional list of window starts to compute; returns amplitude per coarse row */
  function line(cube, xcol, zf0, S, par, win, stat) {
    const no = Math.floor((cube.nz - 1) / SUB) + 1, amp = new Float64Array(no), nb = S.bins.length;
    const recR = S.bins.map(() => new Float64Array(NEL)), recI = S.bins.map(() => new Float64Array(NEL)), yr = new Float64Array(NEL), yi = new Float64Array(NEL);
    const starts = win || (() => { const a = []; for (let w = 0; w < no; w += S.Lw) a.push(w); return a; })();
    for (const w0 of starts) {
      const sp = spectra(cube, w0, S);
      for (let q = 0; q < nb; q++) { const k = S.bins[q]; for (let e = 0; e < NEL; e++) { yr[e] = sp.Xr[e * S.Lp + k]; yi[e] = sp.Xi[e * S.Lp + k]; }
        const Mo = model(xcol, zf0 + w0 * DZ, q, S, par), r = enet(Mo, yr, yi, par); recon(Mo, r.B, recR[q], recI[q]);
        if (stat) { stat.fits++; stat.P += Mo.P; stat.it += r.it; } }
      synth(S, recR, recI, amp, w0);
    }
    return amp;
  }
  return { FULL, FSA, setup, model, enet, recon, spectra, synth, line, erf, colon };
})(M16E);
