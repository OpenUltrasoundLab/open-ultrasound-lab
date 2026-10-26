%% m28_contrast_pulse_sequences.m  --  PI, AM and AM+PI residuals for a bubble and for tissue
% Purpose : transmit the pulse families used in contrast imaging, compute the echo of a
%           lipid-coated bubble (Rayleigh-Plesset/Marmottant, RK4) and of a linear tissue
%           scatterer that sees the same nonlinearly propagated incident pulse, combine the
%           echoes, and measure what survives.
%   PI   : e(+p) + e(-p)          keeps even orders  (2 a2 p^2)
%   AM   : e(p) - 2 e(p/2)        keeps a2 p^2/2 + 3 a3 p^3/4 (odd orders land in the fundamental band)
%   AMPI : e(p) + 2 e(-p/2)       keeps 3 a2 p^2/2 + 3 a3 p^3/4
% Tissue  : incident pulse p_inc = a p1 + a^2 (beta x / (2 rho c^3)) d(p1^2)/dt
%           (second-order plane-wave solution), echo = p_inc delayed by motion 2 v k / (c PRF).
% Change  : P (kPa at the bubble), f, R0, depth x, motion v, half-amplitude error errDb.
% Expect  : at 100 kPa, 2.5 MHz, R0 = 2 um, 40 mm: tissue suppressed by 23.4 dB (PI, limited
%           by the tissue harmonic), 35.4 dB (AM), 25.9 dB (AMPI); bubble retained within 1.5 dB
%           (AM) but 7.5 dB down (PI); 10 mm/s motion at PRF 5 kHz caps AM at about 27 dB;
%           without tissue nonlinearity and motion every sequence cancels tissue exactly.
% Provenance: simulation (this script). Shell values chi = 0.5 N/m, s0 = 0.02 N/m, ks = 5e-9 kg/s,
%           g = 1.07 are illustrative (inside the range fitted to lipid agents, not one agent).
% Runtime : about 6 s.
clear; close all; rng(1);

rho = 1000; cL = 1480; mu = 1e-3; P0 = 101325; sw = 0.072;   % liquid around the bubble
cT = 1540; rhoT = 1000; beta = 4.5;                           % tissue path: c, rho, 1 + B/2A
R0 = 2e-6; g = 1.07; chi = 0.5; ks = 5e-9; s0 = 0.02;         % lipid-coated bubble
f = 2.5e6; N = 3; P = 100e3; x = 0.04; prf = 5000;            % 3-cycle Hann pulse, 40 mm path

Rb = R0 / sqrt(1 + s0/chi); Rr = Rb * sqrt(1 + sw/chi);
sig = @(R) (R >= Rb & R <= Rr) .* chi .* (R.^2/Rb^2 - 1) + (R > Rr) .* sw;
Pg0 = P0 + 2*s0/R0;
Tp = N / f; fs = 100 * f; dt = 1/fs; n = round((Tp + 6/f) / dt) + 1; t = (0:n-1) * dt;
p1  = @(t) (t > 0 & t < Tp) .* 0.5 .* (1 - cos(2*pi*t/Tp)) .* sin(2*pi*f*t);
dp1 = @(t) (p1(t + dt/4).^2 - p1(t - dt/4).^2) / (dt/2);       % d(p1^2)/dt
knl = beta * x / (2 * rhoT * cT^3);                            % s/Pa
fprintf('Tissue 2nd harmonic for a CW wave at this P, f, x: B2/B1 = sigma/2 = %.3f (%.1f dB)\n', ...
  knl * 2*pi*f * P, 20*log10(knl * 2*pi*f * P));

scen = {'default', 0, 0, beta; 'motion 10 mm/s', 0.01, 0, beta; 'half-pulse error +0.1 dB', 0, 0.1, beta; 'linear tissue, no motion', 0, 0, 0};
nfft = 2^nextpow2(4*n); df = fs/nfft;
bandE = @(y, lo, hi) sum(abs(subsref(fft(y, nfft), struct('type', '()', 'subs', {{(ceil(lo/df):floor(hi/df)) + 1}}))).^2);
fprintf('\n%-26s %-5s | tissue: all  f-band  2f-band | bubble: all  f-band  2f-band | CTR gain (dB)\n', 'scenario', 'seq');
for is = 1:size(scen, 1)
  v = scen{is, 2}; errDb = scen{is, 3}; knlS = knl * (scen{is, 4} > 0);
  half = 0.5 * 10^(errDb/20);
  amps = [1, -1, half, -half]; idxp = [0, 1, 1, 1];              % pulse amplitudes and pulse index (motion)
  E_t = zeros(4, n); E_b = zeros(4, n);
  for k = 1:4
    a = amps(k); tau = idxp(k) * 2 * v / (cT * prf) + 1e-7;
    pinc = @(tt) a * P * p1(tt - tau) + a^2 * P^2 * knlS * dp1(tt - tau);
    E_t(k, :) = pinc(t);
    acc = @(tt, R, U) ((Pg0*(R0./R).^(3*g).*(1 - 3*g*U/cL) - P0 - 2*sig(R)./R - 4*mu*U./R - 4*ks*U./R.^2 - pinc(tt))/rho - 1.5*U.^2)./R;
    Rk = R0; Uk = 0; ps = zeros(1, n); ps(1) = rho*(R0^2*acc(0, R0, 0))/0.01;
    for i = 2:n
      tk = t(i-1); h = dt;
      k1R = Uk; k1U = acc(tk, Rk, Uk);
      k2R = Uk + h/2*k1U; k2U = acc(tk + h/2, Rk + h/2*k1R, Uk + h/2*k1U);
      k3R = Uk + h/2*k2U; k3U = acc(tk + h/2, Rk + h/2*k2R, Uk + h/2*k2U);
      k4R = Uk + h*k3U;   k4U = acc(tk + h,   Rk + h*k3R,   Uk + h*k3U);
      Rk = Rk + h/6*(k1R + 2*k2R + 2*k3R + k4R); Uk = Uk + h/6*(k1U + 2*k2U + 2*k3U + k4U);
      ps(i) = rho*(Rk^2*acc(t(i), Rk, Uk) + 2*Rk*Uk^2)/0.01;
    end
    E_b(k, :) = ps;
  end
  W = [1 1 0 0; 1 0 -2 0; 1 0 0 2]; names = {'PI', 'AM', 'AMPI'};
  refT = [bandE(E_t(1,:), 0, fs/2), bandE(E_t(1,:), 0.5*f, 1.5*f)];
  refB = [bandE(E_b(1,:), 0, fs/2), bandE(E_b(1,:), 0.5*f, 1.5*f)];
  for s = 1:3
    rt = W(s, :) * E_t; rb = W(s, :) * E_b;
    T3 = 10*log10([bandE(rt, 0, fs/2)/refT(1), bandE(rt, 0.5*f, 1.5*f)/refT(2), bandE(rt, 1.5*f, 2.5*f)/refT(1)] + 1e-30);
    B3 = 10*log10([bandE(rb, 0, fs/2)/refB(1), bandE(rb, 0.5*f, 1.5*f)/refB(2), bandE(rb, 1.5*f, 2.5*f)/refB(1)]);
    fprintf('%-26s %-5s |       %6.1f %6.1f %7.1f   |       %6.1f %6.1f %7.1f   |  %6.1f\n', scen{is, 1}, names{s}, T3, B3, B3(1) - T3(1));
    if is == 1, keepT(s, :) = rt; keepB(s, :) = rb; end
  end
  if is == 1, E1t = E_t(1, :); E1b = E_b(1, :); end
end
fprintf('(all: residual energy re the single-pulse echo; f-band: 0.5f-1.5f re the single echo in that band; 2f-band re the whole single echo; tissue at -300 dB means exact cancellation)\n');

%% Motion law check: PI residual of a linear scatterer = |1 - exp(-j 2 pi f tau)| = 2 |sin(pi f tau)|
v = 0.01; tau = 2*v/(cT*prf);
fprintf('\nMotion: v = 10 mm/s, PRF 5 kHz -> shift %.2f ns, 2|sin(pi f tau)| = %.1f dB (the f-band floor above)\n', tau*1e9, 20*log10(2*abs(sin(pi*f*tau))));

%% Figure: single echoes and residuals
figure('Position', [50 50 900 520]);
subplot(2, 2, 1); plot(t*1e6, E1t/max(abs(E1t)), 'k', t*1e6, keepT(1,:)/max(abs(E1t)), 'b', t*1e6, keepT(2,:)/max(abs(E1t)), 'r'); grid on;
title('tissue: single echo and residuals'); xlabel('time (\mus)'); ylabel('re single echo'); legend('single', 'PI', 'AM');
subplot(2, 2, 2); plot(t*1e6, E1b/max(abs(E1b)), 'k', t*1e6, keepB(1,:)/max(abs(E1b)), 'b', t*1e6, keepB(2,:)/max(abs(E1b)), 'r'); grid on;
title('bubble: single echo and residuals'); xlabel('time (\mus)'); ylabel('re single echo');
fa = (0:nfft/2-1)*df/1e6; sp = @(y) 20*log10(abs(fft(y, nfft)) + 1e-30);
SpT = sp(E1t); SpB = sp(E1b); mT = max(SpT); mB = max(SpB);
subplot(2, 2, 3); s1 = sp(keepT(1,:)); s2 = sp(keepT(2,:));
plot(fa, SpT(1:nfft/2) - mT, 'k', fa, s1(1:nfft/2) - mT, 'b', fa, s2(1:nfft/2) - mT, 'r'); grid on; xlim([0 10]); ylim([-80 5]);
xlabel('f (MHz)'); ylabel('dB'); title('tissue spectra');
subplot(2, 2, 4); s1 = sp(keepB(1,:)); s2 = sp(keepB(2,:));
plot(fa, SpB(1:nfft/2) - mB, 'k', fa, s1(1:nfft/2) - mB, 'b', fa, s2(1:nfft/2) - mB, 'r'); grid on; xlim([0 10]); ylim([-80 5]);
xlabel('f (MHz)'); ylabel('dB'); title('bubble spectra');
