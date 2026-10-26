# Module 25 beamformer comparison: setup, parameters, results (state Oct 25, 2026)

Single source of truth for the 11-method comparison in Module 25 (Demo 25.6 and related demos). Page: m25.html. Offline scripts: tools/m25-admire-tuning/ (README there). ADMIRE parameter provenance: ADMIRE_NOTES.md.

## Simulated acquisition (all methods see the same channel data)
| Item | Value |
|---|---|
| Array | 64 elements, 0.3 mm pitch, point (omnidirectional) elements |
| Pulse | Gaussian, f0 = 5 MHz, 60 % (−6 dB) fractional bandwidth |
| Sampling | IQ at 20 MHz; image grid 81 × 161 pixels (0.15 × 0.075 mm); DMAS/SLSC use RF on a 4× finer axial grid |
| Transmit | 9 plane waves, −4° to 4° in 1° steps, coherent compounding per channel; transmit grating lobe at 17.6 mm (outside phantom) |
| Receive | rectangular apodization, all 64 elements |
| Phantom | fully developed speckle (80 scatterers/mm², z = 26–34.6 mm); anechoic cyst r = 2.5 mm at (0, 30) mm; bright scatterer (4.2, 30.5) mm; point (−3, 24) mm; pair at x = 2.25 and 2.75 mm, z = 24 mm |
| Displayed seeds | scatterers 7, noise 108 (tuning used held-out seeds 1007/1108) |
| Noise | band-limited complex channel noise (white Gaussian filtered by the pulse spectrum), independent per element/angle/sample, power re mean per-channel speckle power; slider −40 … +30 dB |
| Noise modes | (1) plane-wave shared recordings (physical, default, all panels); (2) line-by-line: independent noise per image line (DAS-only display emulation, or all beamformers via checkbox, ~13 s); (3) per-pixel noise after beamforming (non-physical comparison, DAS only) |
| Image SNR readout | channel 0 dB → DAS speckle-to-noise ≈ +13.9 dB (−20 dB → +33.9; +30 dB → −16.1) |
| Aberration | optional receive phase screen, Gaussian-correlated element delays, 3 mm correlation length, rms slider (ns) |
| Metrics | lateral FWHM and PSL of the point; pair dip; CR (dB, before clipping); CNR at display DR; gCNR (raw and at display DR); speckle SNR; time per pixel |

## Methods, parameters, sources, tuned values
| Method | Default parameters | Source | Tuned (rule below) | Main weakness |
|---|---|---|---|---|
| DAS | uniform, 64 el., 9 angles | classical | default optimal (Hann violates FWHM limit) | sidelobe clutter fills cyst |
| MV | L = 32, K = ±2, Δ = 0.01 | Synnevåg 2007/2009; values ours | L 40, K ±4, Δ 1e−4 (gCNR 0.978 → 0.982) | lone-point width measures SNR; low speckle SNR |
| ESBMV | δ = 0.5 (≤ 4 eigvecs) + MV | Asl & Mahloojifar 2010 | L 48, δ 0.3 (0.988 → 0.991) | breaks wᴴa = 1 next to bright scatterers |
| CF | none | Mallart & Fink 1994 | weight off = DAS (0.935 → 0.958) | dynamic-range stretch, dark holes |
| GCF | M0 = 1 | Li & Li 2003 | default optimal | arbitrary cut-off |
| PCF | γ = 1, σ0 = π/√3 | Camacho 2009 | γ 0.25 (0.914 → 0.957) | clipping to zero at critical γ |
| SLSC | M = 8, kernel ±0.15 mm | Lediju 2011 | M 32, kernel ±0.30 mm (0.991 → 0.963; default violates FWHM limit, 1.87 mm) | coherence, not amplitude; wide point |
| F-DMAS | Gaussian 2f0 band-pass, σ = 1.80 MHz | Matrone 2015 | no free parameter | small gain on broadband data |
| FXPF | p = 4, W = 16 (1.2 mm), μ = 1e−3 | Shin 2016/2018 | p 3, W 32 (0.962 → 0.966) | keeps coherent off-axis clutter |
| Wiener postfilter | per sample | Nilsen & Holm 2010 | no free parameter | weak on noise, none on clutter |
| ADMIRE | α = 0.9, λ = 0.0189·rms, 4-sample rect windows with 50 % overlap (authors' recommendation, Oct 2026; released code: none), 3 bins, ~31 000 predictors (full) / ~900 (reduced) | Byram 2015 + authors' code (VU-BEAM-Lab/ADMIRE) + Schlunk 2021/2023; plane-wave adaptation and calibration ours | default kept (cλ 0.0189 meets FWHM limit with overlap: 0.5263 ≤ 0.5267 mm; held-out gCNR 0.964); before overlap: cλ/4 | cost; exact zeros make CR meaningless |

Tuning rule: maximize gCNR subject to point FWHM ≤ 1.1 × DAS FWHM (1.1 × 0.479 = 0.527 mm) at −20 dB channel noise, no aberration, on the held-out realization; evaluated on the displayed realization. Gains of a few thousandths are within realization-to-realization spread (stated on the page).

## ADMIRE precomputed grid (full model, 50 % STFT overlap since Oct 26)
8 noise levels (−40 … +30 dB, 10 dB steps) × 3 aberration levels (0, 20, 40 ns), shared-recording noise; 10.4 CPU-hours (3.8 without overlap); embedded as `<script id="m25-adm-data">`. Live reduced model ~140 s (Node.js).
| Setting | DAS gCNR · CR | ADMIRE gCNR · CR |
|---|---|---|
| −20 dB, 0 ns | 0.958 · 24.7 dB | 0.988 · 36.0 dB (no overlap: 0.984 · 35.7) |
| 0 dB, 0 ns | 0.815 · 13.6 dB | 0.964 · 27.0 dB (no overlap: 0.942 · 26.5) |
| 0 dB, 40 ns | 0.729 · 9.6 dB | 0.880 · 17.2 dB (no overlap: 0.870 · 17.5) |
| +20 dB, 0 ns | 0.205 · 0.8 dB | 0.213 · 3.2 dB (no overlap: 0.274 · 2.7, inflated by 21 % exact-zero cyst pixels) |

Speckle SNR (background, full model): 1.94 / 1.86 at −20 / 0 dB with overlap vs 1.82 / 1.72 without (DAS 1.86 / 1.88; Rayleigh 1.91). Panel grid: 11 images (Wiener added Oct 26). Comparison with Dei, Tierney & Byram, J Med Imaging 2018 (doi 10.1117/1.JMI.5.2.027001) stated on the page.

## Known limitations / open items
- One simulated phantom, point elements, no reverberation, no attenuation; PICMUS validation not yet done (planned v1.1).
- ADMIRE tuned with the reduced model only; full-model tuning too costly.
- Original ADMIRE 2015/2017 full texts could not be fetched; parameters from authors' code and later papers.
- Line-by-line mode emulates noise statistics only (transmit focusing not simulated); in that mode FWHM/pair-dip columns measure noise.
