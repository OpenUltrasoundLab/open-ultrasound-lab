%% m17_stack_dimensions.m  --  typical array acoustic stack at 3.5 and 7 MHz
% Purpose : compute the dimensions drawn in Figures 17.C-17.E: active-layer
%           thickness (lambda/2 plate), Desilets two-layer impedances and
%           lambda/4 thicknesses, pitch, kerf, element w/t (sub-dicing check),
%           cylindrical lens radius / sag / centre thickness and loss, and the
%           elevation slice thickness and near-field length.
% Model   : rules of thumb only (no wave solver). PZT-5H-like plate
%           (7500 kg/m^3, 4560 m/s), layers c = 2500 m/s, tissue 1540 m/s,
%           RTV lens 1000 m/s with 3 dB/mm at 5 MHz and alpha ~ f^1.5,
%           linear array pitch = lambda. Plate speed is used for t; a slender
%           element resonates lower (m17_fdtd_plate_modes.m).
% Change  : f0 list, H and F per frequency, kerf, cLay, the array pitch rule.
% Expect  : 3.5 MHz: t 0.651 mm, layers 8.95/2.34 MRayl at 179 um, pitch
%           0.440 mm, w/t 0.63, R 32.4 mm, sag 0.560 mm, slice 1.95 mm.
%           7 MHz: t 0.326 mm, layers 89 um, pitch 0.220 mm, w/t 0.61,
%           R 10.8 mm, sag 0.237 mm, slice 0.87 mm. Runtime < 1 s.
clear; close all; rng(1);

%% Parameters
f0s  = [3.5 7];            % MHz
Hs   = [12 4.5];           % elevation height (mm)
Fs   = [60 20];            % elevation focus (mm)
kerf = [0.030 0.020];      % dicing kerf (mm)
vP = 4560; rhoP = 7500; Zp = rhoP*vP/1e6;   % plate speed (m/s), density, MRayl
ZL = 1.5; cLay = 2500;                      % tissue MRayl, layer speed (m/s)
ct = 1540; cl = 1000; te = 0.30;            % tissue / lens speeds (m/s), lens edge (mm)
a5 = 3; nexp = 1.5;                         % lens loss dB/mm at 5 MHz, exponent
Z1 = Zp^(4/7)*ZL^(3/7); Z2 = Zp^(1/7)*ZL^(6/7);   % Desilets two-layer

figure('Position', [80 80 900 420]);
for q = 1:numel(f0s)
  f0 = f0s(q); lam = ct/(f0*1e3);           % mm
  t = vP/(2*f0*1e3);                        % plate lambda/2 (mm)
  d = cLay/(4*f0*1e3);                      % lambda/4 layer (mm)
  pitch = lam; w = pitch - kerf(q);
  R = Fs(q)*(ct/cl - 1); sag = R - sqrt(R^2 - (Hs(q)/2)^2); Tc = te + sag;
  alpha = a5*(f0/5)^nexp;                   % dB/mm at f0
  slice = 0.886*lam*Fs(q)/Hs(q); nf = Hs(q)^2/(4*lam);
  fprintf('--- %.1f MHz (lambda %.3f mm in tissue)\n', f0, lam);
  fprintf('plate t = %.3f mm (Zp %.1f MRayl); ML1 %.2f MRayl, ML2 %.2f MRayl, each %.0f um\n', t, Zp, Z1, Z2, 1e3*d);
  fprintf('pitch %.3f mm, kerf %.0f um, width %.3f mm, w/t = %.2f', pitch, 1e3*kerf(q), w, w/t);
  if w/t > 0.6, fprintf(' -> sub-dice: halves at w/t = %.2f\n', (w - kerf(q))/2/t); else, fprintf(' (ok)\n'); end
  fprintf('lens R = %.1f mm, sag %.3f mm, centre %.3f mm, two-way centre loss %.1f dB\n', R, sag, Tc, 2*alpha*Tc);
  fprintf('elevation F %.0f mm (F# %.2f): slice %.2f mm, near field H^2/(4 lambda) = %.1f mm\n', Fs(q), Fs(q)/Hs(q), slice, nf);
  % stack to scale (azimuth section, 3 elements)
  subplot(1, 2, q); hold on;
  z = 0; lay = {t, d, d, Tc}; col = {[0.85 0.55 0.3], [0.3 0.65 0.4], [0.6 0.85 0.65], [0.75 0.7 0.95]};
  for e = 0:2
    x0 = e*pitch;
    patch(x0 + [0 w w 0], [0 0 t t], col{1});
    patch(x0 + [0 w w 0], t + [0 0 d d], col{2});
    patch(x0 + [0 w w 0], t + d + [0 0 d d], col{3});
  end
  patch([0 3*pitch-kerf(q) 3*pitch-kerf(q) 0], t + 2*d + [0 0 Tc Tc], col{4});
  patch([0 3*pitch-kerf(q) 3*pitch-kerf(q) 0], [-0.2 -0.2 0 0], [0.8 0.8 0.8]);
  axis equal; xlim([-0.05 3*pitch+0.05]); ylim([-0.25 t+2*d+Tc+0.1]);
  xlabel('azimuth (mm)'); ylabel('thickness (mm)');
  title(sprintf('%.1f MHz: plate, 2 layers, lens (to scale)', f0));
end
