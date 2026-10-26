# Module 25 offline computations (ADMIRE precompute + equal tuning)
Reproduces the precomputed data embedded in m25.html.
- `load.js` evaluates the page's own engine (M16E + ADM) in Node, so offline and in-page results use identical code.
- `precompute.js k n` — worker k of n; full-model ADMIRE for 8 noise × 3 aberration settings; resumable; writes out/lines_k.jsonl (3.8 CPU-hours total on 2 cores).
- `merge.js` → admire_grid.json + admire_imgs.json; `build_data.py` embeds them into m25.html (`<script id="m25-adm-data">`).
- `validate_a.js` — point target + off-axis clutter check (log in validate_a.log).
- `calibrate.js` — per-bin wavenumber calibration fit.
- `tune.js`, `gen7c.py` — equal tuning protocol (max gCNR s.t. FWHM ≤ 1.1× DAS, held-out seeds 1007/1108) and trade-off figure data.
- `ADMIRE_NOTES.md` — every parameter with its source ([repo]/[paper]/[ours]).
Paths inside the scripts point to the original work folder; adjust `__dirname`/page path if run elsewhere.

Raw outputs (out/, ~11 MB) are not in the repository; rerun precompute.js to regenerate them. The page uses the merged result embedded in m25.html.
