%% m09_fir_window_kaiser.m
% Module 09 - FIR design without a toolbox: Kaiser window formulas vs equiripple.
%  1) Kaiser: beta(A) and M = (A - 8)/(2.285 dw) for an RF anti-alias low-pass
%     (fs = 40 MHz, pass <= 8 MHz, stop >= 9 MHz, A = 60 dB); measure what it achieves
%  2) Equiripple order estimates (Kaiser 1974, Herrmann-Rabiner-Chan 1973) for
%     dp = 0.01, ds = 0.001, and an equiripple design by Lawson's iteratively
%     reweighted least squares (IRLS), which needs only backslash
% What to change: fp, fst, A, dp, ds, the IRLS orders.
% Expected output: beta = 5.653, M = 145 -> about 59.6 dB; A = 61 -> M = 148, 60.9 dB;
% estimates 101.4 (both); IRLS fails at order 102 and meets dp <= 0.01, ds <= 0.001 from 104.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.

clear; close all; clc;
rng(1);

%% Spec
fs  = 40;            % sampling rate (MHz)
fp  = 8;  fst = 9;   % passband and stopband edges (MHz)
wp  = 2*pi*fp/fs; ws = 2*pi*fst/fs; dw = ws - wp;   % rad/sample
wgrid = linspace(0, pi, 16384)';                      % dense evaluation grid

%% Helpers (anonymous only)
k = 0:40; I0 = @(x) sum(((x(:)/2).^(2*k)) ./ (factorial(k).^2), 2);       % Bessel I0, series
kbeta = @(A) (A > 50) .* 0.1102 .* (A - 8.7) + (A >= 21 & A <= 50) .* (0.5842*max(A-21,0).^0.4 + 0.07886*(A-21));
kord  = @(A, dw) (A - 8) ./ (2.285 * dw);
kwin  = @(M, b) I0(b * sqrt(max(0, 1 - (2*(0:M)'/M - 1).^2))) / I0(b);
ideal = @(M, wc) (wc/pi) * ((0:M)' == M/2) + ((0:M)' ~= M/2) .* sin(wc*((0:M)' - M/2)) ./ (pi*((0:M)' - M/2) + ((0:M)' == M/2));
Hmag  = @(h) abs(exp(-1j * wgrid * (0:numel(h)-1)) * h(:));
ripple_p = @(H) max(abs(H(wgrid <= wp) - 1));
ripple_s = @(H) max(H(wgrid >= ws));

%% 1) Kaiser window design
fprintf('Kaiser window, pass %g MHz, stop %g MHz at fs = %g MHz\n', fp, fst, fs);
designs = {};
for A = [60 61]
    b = kbeta(A); Mf = kord(A, dw); M = ceil(Mf);
    h = ideal(M, (wp + ws)/2) .* kwin(M, b);
    H = Hmag(h); dp_ = ripple_p(H); ds_ = ripple_s(H);
    fprintf('  A = %d dB: beta = %.3f, M = %.1f -> %d (%d taps); measured stopband %.2f dB, passband dev %.4f\n', ...
        A, b, Mf, M, M+1, 20*log10(ds_), dp_);
    designs{end+1} = h; %#ok<SAGROW>
end
hK = designs{2}; MK = numel(hK) - 1;
fprintf('  delay M/2 = %d samples = %.3f us = %.2f mm (c/2 = 770 m/s)\n', MK/2, MK/2/fs, MK/2/fs*0.770);

%% 2) Equiripple order estimates
dp = 0.01; ds = 0.001; df = (fst - fp)/fs;
Nk = (-20*log10(sqrt(dp*ds)) - 13) / (14.6*df);
a6 = [0.005309 0.07114 -0.4761 0.00266 0.5941 0.4278]; lp = log10(dp); ls = log10(ds);
Dinf = ls*(a6(1)*lp^2 + a6(2)*lp + a6(3)) - (a6(4)*lp^2 + a6(5)*lp + a6(6));
fK = 11.01217 + 0.51244*(lp - ls); Nh = Dinf/df - fK*df;
fprintf('\nEquiripple estimates (dp = %.3f, ds = %.4f, df = %.4f): Kaiser %.1f, Herrmann %.1f\n', dp, ds, df, Nk, Nh);
fprintf('  A/22 rule: %.0f taps for A = %d dB\n', -20*log10(ds)/(22*df), round(-20*log10(ds)));

%% 3) Equiripple by Lawson IRLS (type I linear phase, M even)
% A(w) = sum_k a_k cos(k w); weighted LS on a grid of the two bands, weights updated
% as u <- u .* |E| so the error envelope flattens toward the minimax solution.
Wst = dp/ds;                                          % stopband weight
fprintf('\nIRLS equiripple, weight %.0f in the stopband:\n', Wst);
hE = [];
for M = [102 104 106 108]
    L = M/2; ng = 16*M;
    gp = linspace(0, wp, round(ng*wp/(pi - dw)))'; gs = linspace(ws, pi, round(ng*(pi - ws)/(pi - dw)))';
    g = [gp; gs]; d = [ones(size(gp)); zeros(size(gs))]; W = [ones(size(gp)); Wst*ones(size(gs))];
    C = cos(g * (0:L)); u = W;
    for it = 1:80
        a = (C .* u) \ (d .* u);                     % weighted least squares
        e = abs(W .* (C*a - d));
        u = u .* sqrt(e + 1e-14); u = u / sqrt(mean(u.^2));
    end
    h = [flipud(a(2:end))/2; a(1); a(2:end)/2];
    H = Hmag(h); dp_ = ripple_p(H); ds_ = ripple_s(H);
    ok = dp_ <= dp && ds_ <= ds;
    fprintf('  M = %d: passband dev %.5f, stopband %.6f (%.2f dB) -> %s\n', M, dp_, ds_, 20*log10(ds_), char('fails' .* ~ok + 'meets' .* ok));
    if M == 106, hE = h; end
end
fprintf('  (true Remez minimum for this spec: order 104; equiripple with dp = ds = 0.001 needs 132)\n');

%% Figures
f = wgrid/(2*pi)*fs;
figure('Name', 'Kaiser vs equiripple');
subplot(2,1,1); plot(f, 20*log10(Hmag(hK) + 1e-12), 'b', 'LineWidth', 1); hold on;
plot([fst fs/2], [-60 -60], 'r--', [fp fp], [-120 5], 'r:', [fst fst], [-120 5], 'r:');
ylim([-110 5]); grid on; xlabel('frequency (MHz)'); ylabel('|H| (dB)'); title(sprintf('Kaiser window, M = %d, beta = %.2f', MK, kbeta(61)));
subplot(2,1,2); plot(f, 20*log10(Hmag(hE) + 1e-12), 'Color', [0.8 0.2 0.4], 'LineWidth', 1); hold on;
plot([fst fs/2], [-60 -60], 'r--', [fp fp], [-120 5], 'r:', [fst fst], [-120 5], 'r:');
ylim([-110 5]); grid on; xlabel('frequency (MHz)'); ylabel('|H| (dB)'); title('Equiripple (IRLS), M = 106, dp = 0.01, ds = 0.001');
