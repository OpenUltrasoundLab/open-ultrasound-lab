# ADMIRE in module 25 — implementation notes

Status log at the bottom. Every parameter has a source tag: [repo] = authors' code, [paper:X] = a fetched paper, [ours] = our choice (with reason).

## 1. Sources actually read (Oct 24)
- [repo] VU-BEAM-Lab/ADMIRE, GitHub (Apache-2.0), cloned Oct 24 2026. Authors' own MATLAB/C/CUDA code (Khan, Dei, Schlunk, Byram).
  Files read: README.md, ADMIRE_models_generation_main.m, generate_models_for_stft_window.m, generate_model_for_frequency.m,
  generate_modeled_signal_for_predictor.m, wavenumber_calibration_value_lookup.m, apply_ADMIRE_models_CPU.m, ccd_double_precision.c.
  README states the code is "meant for walked aperture scans with focused transmits" and cites Byram et al. 2015 for the method.
- [paper:Schlunk2021] Schlunk, Dei, Byram, "Iterative model-based beamforming for high dynamic range applications," IEEE TUFFC 68(3):482–493, 2021, doi 10.1109/TUFFC.2020.3012165 (NSF PAR full text). α = 0.9; λ = (0.0189/2)·y_rms (LDF) or (0.0189/10)·y_rms (HDF); 90 % STFT overlap; in vivo: 8 plane waves −9…9° compounded into a synthetically focused data set, then ADMIRE.
- [paper:SchlunkIUS2021] Schlunk & Byram, "On model space sampling in ADMIRE for image quality and computational efficiency," IEEE IUS 2021, doi 10.1109/IUS52206.2021.9593649 (NSF PAR). Base model 27 574 predictors; λ = 0.00189·y_rms, α = 0.9. Cyst (no added clutter): DAS CR −32.2 ± 1.0 dB, CNR 6.1 ± 0.6; ADMIRE base model CR −44.3 ± 1.3 dB, CNR 5.2 ± 0.6. With 0 dB SCR clutter: DAS CR −12.1 dB, CNR 3.5; ADMIRE (Base Rand) CNR 4.3.
- [paper:Schlunk2023] Schlunk, Hsi, Byram, Ultrasonics 131:106952, 2023, doi 10.1016/j.ultras.2023.106952 (NSF PAR). ADMIRE applied to curvilinear plane-wave synthetic focusing (64 angles); α = 0.9, λ = 0.00189·y_rms "based on previous work". Their Table 1 model sampling (ROI 0.179 res_lat × 0.715 res_ax, clutter 3.577 res_lat × 7.154 res_ax) is sparser than the repo default.
- [paper:Khan2020] Khan et al., "Real-time, simultaneous DAS, ADMIRE and SLSC imaging using GPU-based processing," IEEE IUS 2020 (NSF PAR): padded STFT window 12 samples, 0 % overlap, 3 frequencies, α = 0.9, c_λ = 0.0189, tolerance 0.1 with ICA, focused walked aperture (L7-4).
- [paper:Schlunk2019iADMIRE] "Iterative ADMIRE for high dynamic range B-mode" (NSF PAR, IUS): α = 0.9; CR dynamic range DAS 50.4 dB, ADMIRE 49.3 dB, iADMIRE 69.6 dB.
- Crossref metadata verified: Byram, Dei, Tierney, Dumont, IEEE TUFFC 62(11):1913–1927, 2015, doi 10.1109/TUFFC.2015.007004; Dei & Byram, IEEE TUFFC 64(10):1450–1464, 2017, doi 10.1109/TUFFC.2017.2729944.
- NOT READ (blocked: PMC/Europe PMC/IEEE/SPIE returned captcha, 403 or rate limits on every attempt): full texts of Byram 2015 (PMC4778405) and Dei & Byram 2017 (PMC5656016), Dei & Byram SPIE 2017 "ADMIRE with plane wave synthesis" (doi 10.1117/12.2255526), Dei & Schlunk & Byram TUFFC 2019 (doi 10.1109/TUFFC.2019.2924824). Consequence: parameters below come from the authors' released code (which implements the 2015 model) and from the later papers above; we do not quote numbers from the 2015/2017 papers.

## 2. The algorithm as released [repo]
1. Delay the channel data (dynamic receive focusing for the beam).
2. STFT along depth per channel: rectangular window, length L_w = ceil(pulse_multiplier · FWHM_pulse · f_s), FWHM_pulse = sqrt(8 ln2)·σ_t, σ_t = 1/(2π σ_f), σ_f = BW·f0/sqrt(8 ln2); pulse_multiplier = 1; zero-padded to 2 L_w; overlap 0 ("must set to 0 in this code release"); window function @rectwin (must).
3. Frequencies fitted: FFT bins with f0 − 0.6·BW·f0 ≤ f ≤ f0 + 0.6·BW·f0 (positive half only). Others set to 0 in the reconstruction; negative bins = conj (real RF).
4. Model per window (centre depth z_c = mean depth of the window) and per frequency: columns = modelled aperture signals of scatterers at (x, z) with a distance offset d (extra path, models multipath/arrival time):
   modelled signal = A_FT · A_es · exp(i k zd), zd = 2 z' where z' solves sqrt(x_e^2 + z'^2) + z' = sqrt((x_e − x)^2 + (z − cal)^2) + d + d_shift + 2 z_c − (z − cal);
   A_FT = sqrt(normalized integral of a Gaussian power envelope exp(−((u − zd)/s_t)^2), s_t = c σ_t, over [max(zd − h, 2 z_first), min(zd + h, 2 z_last)]), h = L_w/f_s · c/2; zero if outside;
   A_es = sinc(π p sinθ/λ)·cosθ (element factor, element width = pitch).
   - resolution scales: res_lat from the FWHM of the FFT (4096 pt) of |on-axis model| → (fwhm bins)·(1/p)/4096·λ·z_c; res_ax = 2 res_lat.
   - ROI model: x ∈ [−0.5 : 0.0179 : 0.5]·res_lat, z ∈ z_c + [−0.5 : 0.143 : 0.5]·res_ax, d ∈ [−0.8 mm : 0.0485 λ : 0.4 mm]; keep (x/res_lat)^2 + ((z − z_c)/res_ax)^2 ≤ 1.
   - Outer model: x ∈ [−L : 1.4228 res_lat : L], L = p·N/2 + 1 mm; z ∈ [0 : 0.7114 res_ax : 1.05 z_c] (default) or [z_c − 0.5 res_ax : 0.7114 res_ax : z_c + 0.5 res_ax] (README's documented off-axis-only alternative, "tested"); d ∈ [−8 mm : 0.1211 λ : 3.2 mm]; keep outside the ellipse with semi-axes res_lat + 0.5 mm, res_ax + 0.5 mm.
   - Predictors whose signal is all zero are dropped. Each complex column normalized to unit L2 norm. Optional ICA (we do not use it; 2015 "full model").
   - Calibration constants: cal_shift 7.75 µm, distance_offset_shift 38.5 µm, wavenumber calibration factors per bin (e.g. [1.29 0.995 0.85 0.835] at 5 MHz) — empirical for their probes.
5. Real-ification: X_r = [[Re X; Im X], [−Im X; Re X]], y = [Re y; Im y] → each predictor gets a complex coefficient.
6. Elastic net (cyclic coordinate descent, C code; glmnet originally):
   β̂ = argmin (1/2N)||y − Xβ||^2 + λ(α||β||_1 + (1−α)||β||^2/2), N = 2·channels.
   y standardized by its (population) std before the fit, λ divided by the same std, β rescaled after.
   λ = c_λ · sqrt(mean(y^2)), c_λ = 0.0189, α = 0.9; update β_j = S(p_j, λα)/(1/N + λ(1−α)), p_j = x_j^T r_(−j)/N.
   Stop when max_j (Δβ_j)^2/N < tol; tol = 0.1 with ICA, "1e-4 or 1e-5" without ICA (README); max 100 000 iterations.
7. Reconstruction: y_ROI = X_ROI β_ROI (both real-ified halves of the ROI predictors), back to complex, placed in the selected bins (conj in mirrored bins), IFFT, first L_w samples → decluttered channel data; then sum over elements, envelope. Depths without ADMIRE use DAS.

## 3. Our implementation and adaptations
Code: `const ADM = (function (E) {...})(M16E);` in m25.html right after the engine (no DOM); the page and the offline script (Node, eval of the page's own M16E + ADM + GRID/PROF code) run the same text.

| Item | Value used | Source |
|---|---|---|
| Input data | delayed, phase-rotated, 9-angle compounded aperture data of Demo 25.6 (`M16E.column`), all 64 elements, coarse image rows (Δz = 0.075 mm) | ours: plane-wave adaptation = compound first, then ADMIRE on the synthetically focused aperture data, as in [paper:Schlunk2021] (8 PW, in vivo) and [paper:Schlunk2023] (curvilinear PW synthetic focusing). Release code is for focused walked apertures. |
| Axial sampling for the STFT | complex (analytic) samples on the image grid, f_s,eq = c/(2Δz) = 10.27 MHz | ours: same grid as FXPF; gives the release's typical 3 bins around f0 |
| STFT window | rectangular, L_w = ceil(1·FWHM·f_s) = ceil(0.294 µs · 10.27 MHz) = 4 samples (0.30 mm), zero-padded to 8, 50 % overlap (hop 2; overlapping windows averaged per row) since Oct 26, 2026 [ours: ADMIRE authors' recommendation, personal communication; released code has no overlap; reduced-model test: speckle SNR 1.77→1.91 at −20 dB] | [repo] rule (pulse_multiplier 1, rectwin, overlap 0) |
| Bins fitted | f0 ± 0.6·BW·f0 = 3.2–6.8 MHz → 3.85, 5.13, 6.42 MHz | [repo] rule; ours: no positive-half restriction because the data are analytic (no mirrored negative band) |
| Wavenumber calibration | k_eff = k · [1.10, 0.95, 0.88] per bin | ours, fitted like the release's per-probe factors: maximize mean model–data correlation for 42 isolated point scatterers (2 lines, 7 lateral offsets 0–5 mm, 3 depth offsets). Mean corr 0.882→0.898, 0.905→0.913, 0.824→0.900; in-ROI corr ≥ 0.999 after calibration. Same trend as the release's 5 MHz factors [1.29, 0.995, 0.85, 0.835] (bins above f0 < 1, below > 1). |
| Element factor A_es | 1 | ours: the demo simulates point (omnidirectional) elements |
| Calibration offsets cal_shift, distance_offset_shift | 0 | ours: release values (7.75 µm, 38.5 µm) are empirical for its probes; our delays are exact (model–data corr 0.9998–1.0000 for an on-axis point) |
| Phase sign | exp(−j k zd) | ours, matches the e^{−j2πft} forward FFT and our IQ convention; verified: corr 0.991 vs 0.650 with the opposite sign for points off the line |
| Model space ROI | x ∈ [−0.5:0.0179:0.5]·res_lat, z ∈ z_c+[−0.5:0.143:0.5]·res_ax, d ∈ [−0.8 mm : 0.0485λ : 0.4 mm], ellipse | [repo] defaults |
| Model space outer | x ∈ [−L:1.4228·res_lat:L], L = 64·0.3/2 + 1 = 10.6 mm (centred on the image line); z ∈ z_c + [−0.5:0.7114:0.5]·res_ax; d ∈ [−8 mm : 0.1211λ : 3.2 mm]; outside ellipse (res+0.5 mm) | [repo] defaults except the depth range = README's documented off-axis-only alternative [1, −0.5, 0.7114, 1, 0.5, 0] (ours: the simulation has no reverberation, and the default 0…1.05 z_c depth range gives 150 000–460 000 predictors per fit, measured, infeasible here) |
| res_lat, res_ax | from FFT(4096) of the on-axis model magnitude, res_ax = 2 res_lat (0.45–0.86 mm lateral at our depths/bins) | [repo] |
| Predictors per fit (full) | ≈ 31 000 complex (ROI ≈ 19 600–32 200 by bin, outer ≈ 3 300–13 000) | measured |
| Aperture growth | off (all 64 elements, F/1.3–F/1.9) | ours: every panel of Demo 25.6 uses the same 64-element receive aperture |
| ICA | off ("full model" in the authors' terminology) | [repo] option |
| Elastic net | α = 0.9, λ = c_λ·rms(y), c_λ = 0.0189, standardized y, unit-norm columns, real-ified complex | [repo] |
| Solver | cyclic coordinate descent, zero start, sweeps over the active set between full sweeps; stop when a full sweep changes no coefficient by more than tol | ours (same fixed point as the release's plain cyclic sweeps; checked) |
| Tolerance | 1e-5 on max (Δβ)²/N | [repo] README: "0.0001 or 0.00001" without ICA. Check vs converged (1e-9) solution, 12 fits: ROI reconstruction within 0.5–3.3 % (active set, 1e-5) vs 1.9–12 % for plain cyclic at 1e-4 |
| Warm starts | not used (single λ per fit; warm starts across noise levels changed the stopped solution by 11–16 % at tol 1e-4 → would make results order-dependent) | ours |
| Reconstruction | X_ROI β_ROI per bin, other bins 0, 8-point IFFT, first 4 samples, sum over elements, |·|/64 | [repo] |
| Depth coverage | lines computed over 164 rows (41 full windows), 161 kept, so no DAS fill-in rows | ours |
| Profile lines (FWHM, PSL, pair dip) | one window = rows 14–17 of the 33-row profile line (target row 16 inside it) | ours |
| Reduced model (live) | ROI x/z/d sampling 8×/4×/8× coarser ([−0.5:0.1432:0.5], [−0.5:0.572:0.5], 0.388λ), outer d 8× coarser (0.9688λ); same bins, α, λ, tol; ≈ 900 predictors per fit | ours; 4 test lines at −20 dB: 2.1 dB rms difference from full (pixels within 40 dB of max); full: 43/53 cyst pixels exactly 0, reduced 53/53 |

## 4. Validation results
### (a) Point target + off-axis point clutter (work/validate_a.js, page code, line x = 0)
Scene: unit point at (0, 30) mm; five points of amplitude 3 at x = 1.5, 2.5, −2, −3.5, 4.5 mm, z = 29.95–30.1 mm; complex pulse-band noise at −40 dB re a unit echo's peak power.
| | full model | reduced model |
|---|---|---|
| target alone: ADMIRE/DAS at the target pixel | −5.0 dB (shrinkage + 3-bin band; irrelevant after normalization) | −5.3 dB |
| clutter only, at the target pixel (re each method's target-alone value) | DAS −33.6 dB, ADMIRE exactly 0 | DAS −33.6 dB, ADMIRE exactly 0 |
| clutter only, rms over rows > 0.45 mm from the target depth | DAS −40.1 dB, ADMIRE −64.0 dB | DAS −40.1, ADMIRE −84.4 dB |
| target + clutter at the target pixel (re own target-alone) | DAS 0.0 dB, ADMIRE −0.7 dB | DAS 0.0, ADMIRE −0.7 dB |
| lateral FWHM of the lone target (z = 30.025 mm) | DAS 0.582 mm, ADMIRE 0.628 mm | ADMIRE 0.671 mm |
Verdict: keeps the point (−0.7 dB with clutter present), removes the off-axis clutter at the target pixel completely and by 24 dB elsewhere. The mainlobe is not narrowed (ROI = ±0.5 res_lat; no outer predictors within res_lat + 0.5 mm laterally, so energy from 0.3–1.06 mm off the line can only be fitted by ROI columns and stays).

## 5. Page changes (m25.html)
- 25.5: new h3 "Model-based suppression: ADMIRE" (id m25-admire) replacing the old one-paragraph callout: eqs (25.13) model, (25.14) elastic net + CD step + reconstruction, parameter/source table, Demo 25.5b (one window: aperture fit, selected predictors, depth profile; λ and model-size controls), Demo 25.5c (8 × 3 small multiples, click opens setting in 25.6; gCNR, zero-fraction, speckle-SNR curves).
- Equation renumbering: old (25.13) → (25.15), old (25.14) → (25.16), ids m25-eq13/14 → m25-eq15/16; all eqrefs in m25 updated. No other module referenced them (grep); AUDIT_m25.md keeps the old numbers (historical).
- 25.6: 10th canvas m25-g-admire (grid now 5 columns ≥ 760 px, 2 columns below), readout m25-g-adm-o, button m25-g-adm-run (live reduced), ADMIRE row in the metrics table, h2 "ten beamformers and ADMIRE".
- Engine: `ADM` block after M16E; precomputed data embedded as `<script type="application/json" id="m25-adm-data">` (metrics JSON + PNG sprite, 8 noise columns × 3 aberration rows of 81 × 161 tiles). build_site.py inlines module pages, so no separate asset is needed.
- Rule card (ADMIRE), cheat-sheet row, interview I14, self-check Q9, refs (Byram 2015 with doi, Dei & Byram 2017, Khan 2021 + repo, Schlunk 2021). MATLAB: matlab/m25_admire.m.

## Status log
- Oct 24: literature + repo read; notes started.
- Oct 24 21:57: ADM embedded in m25.html (after M16E). Offline precompute launched: work/precompute.js 2 workers (scratchpad m25adm/work/out/lines_{0,1}.jsonl, resumable: rerun `node precompute.js k 2`). 497 lines (81 grid + 215 profile + 201 pair) x 24 settings.
- Oct 25 (7c): equal tuning: ADMIRE c_λ searched with the REDUCED model only ({0.0189/4, /2, 1, ·2}, α = 0.9): held-out pick 0.0189/4 (default FWHM 0.5271 mm fails the 0.5267 mm ceiling); displayed reduced gCNR 0.982 → 0.976. Pair-dip n/a for precomputed settings +30/20, +20/40, +30/40. Live reduced run measured 76 s (default c_λ) / 116 s (tuned) in headless Chromium.


## Change log
- 2026-10-26: 50 % STFT overlap adopted (ADM.FULL.overlap = 0.5) after feedback from the ADMIRE authors. Full grid recomputed (10.4 CPU-h). Re-tuning keeps c_λ = 0.0189. Old no-overlap outputs: out_noovl_2026-10-25.
