%% m28_ulm_concentration.m  --  ULM: isolation probability, optimal density, capillary waiting time
% Purpose : bubbles in the image plane are a Poisson process of areal density n (per mm^2).
%           A bubble is usable when no other bubble lies inside the ellipse
%           (dx/(k Fx))^2 + (dz/(k Fz))^2 < 1 (Fx, Fz = PSF FWHM). Then
%             P_iso = exp(-n A),  A = pi k^2 Fx Fz,  isolated per mm^2 = n exp(-n A),
%           maximal at n* = 1/A with P_iso = 1/e. Monte Carlo (periodic box) checks this.
%           Second constraint: a capillary of diameter d and speed v receives bubbles at the
%           rate C Q = C (pi d^2/4) v, so its first detection waits about 1/(C Q).
% Change  : f, PSF widths (in wavelengths), k, blood concentration C, slice thickness, vessel table.
% Expect  : Monte Carlo P_iso within about 0.01 of exp(-nA); at 5 MHz (Fx = 1.0, Fz = 0.75 lambda, k = 1)
%           A = 0.224 mm^2 and n* = 4.5 per mm^2; at C = 140 per uL an 8 um capillary at
%           1 mm/s waits about 140 s between bubbles.
% Provenance: Poisson Monte Carlo checked against the exact exp(-nA) law.
% Runtime : about 6 s in Octave.
clear; close all; rng(1);

c = 1540; f = 5e6; lam = c / f * 1e3;           % mm
Fx = 1.0 * lam; Fz = 0.75 * lam; k = 1;          % PSF FWHM (mm) and exclusion factor
A = pi * k^2 * Fx * Fz;                          % exclusion area (mm^2)
fprintf('lambda = %.3f mm, Fx = %.3f mm, Fz = %.3f mm, exclusion area A = %.4f mm^2, n* = 1/A = %.2f per mm^2\n', lam, Fx, Fz, A, 1/A);

%% Monte Carlo in a periodic 10 x 10 mm box
L = 10; ns = [0.5 1 2 4.47 8 16];                % per mm^2
fprintf('\n  n (/mm^2) |  nA   | P_iso MC | exp(-nA) | isolated per mm^2 (MC)\n');
nTrial = 60; Pmc = zeros(size(ns)); Iso = zeros(size(ns));
for in = 1:numel(ns)
  cnt = 0; tot = 0;
  for tr = 1:nTrial
    N = round(ns(in) * L^2);          % fixed count: P_iso = (1 - A/L^2)^(N-1), which is exp(-nA) to < 0.1 %
    xy = rand(N, 2) * L;
    dx = abs(xy(:, 1) - xy(:, 1)'); dx = min(dx, L - dx);    % periodic distances
    dz = abs(xy(:, 2) - xy(:, 2)'); dz = min(dz, L - dz);
    inside = (dx / (k*Fx)).^2 + (dz / (k*Fz)).^2 < 1;
    inside(1:N+1:end) = false;
    cnt = cnt + sum(~any(inside, 2)); tot = tot + N;
  end
  Pmc(in) = cnt / tot; Iso(in) = cnt / (nTrial * L^2);
  fprintf('  %8.2f  | %5.2f |  %.4f  |  %.4f  |  %.3f\n', ns(in), ns(in)*A, Pmc(in), exp(-ns(in)*A), Iso(in));
end
fprintf('Theoretical maximum: %.3f isolated per mm^2 at n* (= 1/(e A))\n', 1/(exp(1)*A));

%% From blood concentration to image density, and the capillary waiting time
C = 140;               % bubbles per uL (= per mm^3); about 2.4 mL x 3e8/mL diluted in 5 L
slice = 0.5;           % mm, elevation slice thickness
vf = 0.05;             % vascular volume fraction
n = C * slice * vf;    % per mm^2 of image
fprintf('\nC = %d per uL, slice %.1f mm, vascular fraction %.0f %%: n = %.2f per mm^2, nA = %.2f, P_iso = %.2f\n', C, slice, 100*vf, n, n*A, exp(-n*A));
ves = [8 1; 20 3; 50 10; 200 50];                   % diameter (um), speed (mm/s)
fprintf('\n vessel d (um) | v (mm/s) | Q (nL/s) | bubbles per s at C | mean wait (s)\n');
for iv = 1:size(ves, 1)
  d = ves(iv, 1) * 1e-3; v = ves(iv, 2);            % mm, mm/s
  Q = pi * d^2 / 4 * v;                             % mm^3/s
  fprintf('   %6.0f      |  %5.1f   | %8.4f |     %9.4g     | %9.3g\n', ves(iv, 1), v, Q * 1e3, C * Q, 1 / (C * Q));
end

%% Doppler shift of a slow bubble: why SVD cannot separate capillary flow from tissue motion
for v = [0.5 1 10] * 1e-3
  fprintf('Axial velocity %.1f mm/s at 5 MHz: Doppler %.1f Hz; displacement per frame at 1 kHz = %.1f um\n', v*1e3, 2*v*f/c, v/1000*1e6);
end

%% Figure
nn = logspace(-1, 2, 200);
figure('Position', [50 50 900 330]);
subplot(1, 2, 1); semilogx(nn, exp(-nn*A), 'b', ns, Pmc, 'ko'); hold on;
semilogx(nn, nn .* exp(-nn*A) * A * exp(1), 'r'); grid on;
xlabel('bubbles per mm^2 of image'); ylabel('fraction');
legend('P_{iso} theory', 'P_{iso} Monte Carlo', 'isolated per area / max', 'Location', 'west'); title('5 MHz, k = 1');
subplot(1, 2, 2); Cs = logspace(0, 4, 200); d = 8e-3; v = 1; Q = pi*d^2/4*v;
loglog(Cs, 1 ./ (Cs * Q), 'b'); grid on; xlabel('blood concentration (bubbles per \muL)'); ylabel('mean wait (s)');
title('8 \mum capillary at 1 mm/s');
