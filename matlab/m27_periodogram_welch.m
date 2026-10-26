%% m27_periodogram_welch.m  -  periodogram bias/variance, chi-square statistics, Welch averaging
% Purpose : speckle-like RF (white Gaussian scatterers convolved with a 5 MHz,
%           60 % bandwidth Gaussian pulse, fs = 40 MHz) has a known power
%           spectral density S(f). Show that
%           (1) the periodogram P(f) is unbiased for smooth S but its relative
%               standard deviation stays 1 however long the record (not
%               consistent); P/S is exponential (chi-square with 2 dof / 2);
%           (2) on a dB scale one bin has std 5.57 dB and sits 2.51 dB low;
%           (3) Welch averaging (Hann, 50 % overlap) cuts the variance to
%               (1 + 2 sum (1-m/Nb) rho(m)^2)/Nb, rho = taper overlap correlation,
%               while the bin width grows from fs/N to fs/M;
%           (4) chi-square confidence intervals with nu = 2 K_eff cover the truth
%               95 % of the time.
% Change  : N (record length), M (Welch segment length), number of trials.
% Expect  : std(P/S) about 1.00 for N = 64, 256, 1024; dB mean -2.5, std 5.6;
%           Welch N = 1024, M = 128: Nb = 15, variance 0.070 predicted, 0.074 measured;
%           CI coverage about 0.95 for both.
clear; close all; rng(1);

%% Parameters
fs = 40;  f0 = 5;  bw = 0.6;                 % MHz; -6 dB fractional bandwidth (amplitude)
sigf = bw*f0/(2*sqrt(2*log(2)));             % std of the amplitude spectrum (MHz)
Nfft = 4096;                                 % length of the generated stationary record
fF = (0:Nfft-1)'/Nfft*fs;  fF(fF >= fs/2) = fF(fF >= fs/2) - fs;
H = exp(-(abs(fF) - f0).^2/(2*sigf^2));      % amplitude transfer function (two-sided, real, even)
Strue = @(f) exp(-(abs(f) - f0).^2/sigf^2);  % PSD of x for unit-variance white input (per sample, fs = 1 units)
% x = ifft(fft(white).*H) has PSD |H|^2 in "per bin" units; P = |X|^2/N estimates it.
makeX = @(n) real(ifft(fft(randn(Nfft, n)).*repmat(H, 1, n)));
inband = @(f) abs(f - f0) <= 1.5;            % +/-1.5 MHz around f0 (inside the -14 dB band)

%% 1. Periodogram: variance does not fall with N
ntr = 200;
fprintf('Periodogram (rectangular window), %d realizations, bins within f0 +/- 1.5 MHz:\n', ntr);
fprintf('    N    bin width (MHz)   mean(P/S)   std(P/S)   mean dB   std dB\n');
for N = [64 256 1024]
  X = makeX(ntr);  x = X(1:N, :);
  P = abs(fft(x)).^2/N;                      % periodogram, per-bin units
  f = (0:N-1)'/N*fs;  k = inband(f) & f < fs/2;
  R = P(k, :)./repmat(Strue(f(k)), 1, ntr);  % should be Exp(1)
  RdB = 10*log10(R(:));
  fprintf('%5d   %10.3f       %8.3f   %8.3f   %7.2f   %6.2f\n', N, fs/N, mean(R(:)), std(R(:)), mean(RdB), std(RdB));
  if N == 1024, Rkeep = R(:); end
end
fprintf('Theory: mean 1, std 1; dB mean = -gamma*10/ln10 = %.2f dB, dB std = (pi/sqrt 6)*10/ln10 = %.2f dB\n', ...
        -0.5772156649*10/log(10), pi/sqrt(6)*10/log(10));

%% 2. Welch: Hann segments of length M, 50 % overlap
N = 1024;  M = 128;  S = M/2;  starts = 1:S:(N - M + 1);  Nb = numel(starts);
w = 0.5 - 0.5*cos(2*pi*((0:M-1)' + 0.5)/M);  U = sum(w.^2);
rho = zeros(Nb-1, 1);
for m = 1:Nb-1
  sh = m*S;  if sh < M, rho(m) = sum(w(1:M-sh).*w(1+sh:M))/U; end
end
varTheory = (1 + 2*sum((1 - (1:Nb-1)'/Nb).*rho.^2))/Nb;
X = makeX(ntr);  fW = (0:M-1)'/M*fs;  kW = inband(fW) & fW < fs/2;
Pw = zeros(M, ntr);
for b = 1:Nb
  seg = X(starts(b):starts(b)+M-1, :).*repmat(w, 1, ntr);
  Pw = Pw + abs(fft(seg)).^2/U;
end
Pw = Pw/Nb;
% expected value = S convolved with the Hann spectral window; S is smooth, so compare with S directly
Rw = Pw(kW, :)./repmat(Strue(fW(kW)), 1, ntr);
fprintf('\nWelch N = %d, M = %d (Hann, 50%% overlap): Nb = %d segments, rho(M/2) = %.4f\n', N, M, Nb, rho(1));
fprintf('  bin width fs/M = %.3f MHz (vs %.3f MHz raw)\n', fs/M, fs/N);
fprintf('  variance of P/S: measured %.4f, theory %.4f, K_eff = %.1f, 1/Nb = %.4f\n', var(Rw(:)), varTheory, 1/varTheory, 1/Nb);
fprintf('  mean of P/S = %.3f (window-smoothing bias is small for this smooth spectrum)\n', mean(Rw(:)));

%% 3. Chi-square confidence intervals and their coverage
% 95 % quantiles of chi2 with nu dof via the Wilson-Hilferty approximation (base MATLAB has no chi2inv)
z975 = 1.959964;
chi2q = @(nu, z) nu*(1 - 2/(9*nu) + z*sqrt(2/(9*nu))).^3;
% exact for nu = 2: chi2_2 quantile = -2 ln(1-p)
lo2 = -2*log(1 - 0.025);  hi2 = -2*log(1 - 0.975);
cov2 = mean(2*Rkeep >= lo2 & 2*Rkeep <= hi2);
fprintf('\nRaw periodogram (nu = 2): 95%% interval for S is [P*2/%.3f, P*2/%.4f] = [%.2f dB, +%.2f dB] around P\n', ...
        hi2, lo2, 10*log10(2/hi2), 10*log10(2/lo2));
fprintf('  span %.1f dB; empirical coverage %.3f\n', 10*log10(hi2/lo2), cov2);
nu = 2/varTheory;  lo = chi2q(nu, -z975);  hi = chi2q(nu, z975);
covW = mean(nu*Rw(:) >= lo & nu*Rw(:) <= hi);
fprintf('Welch (nu = 2 K_eff = %.1f): interval [%.2f dB, +%.2f dB], span %.1f dB; coverage %.3f\n', ...
        nu, 10*log10(nu/hi), 10*log10(nu/lo), 10*log10(hi/lo), covW);

%% Figures
x = X(:, 1);  P1 = abs(fft(x(1:N))).^2/N;  f1 = (0:N-1)'/N*fs;
Pw1 = zeros(M, 1);
for b = 1:Nb, Pw1 = Pw1 + abs(fft(x(starts(b):starts(b)+M-1).*w)).^2/U; end
Pw1 = Pw1/Nb;
figure('Position', [100 100 900 360]);
subplot(1, 2, 1);
h1 = f1 < fs/2;  h2 = fW < fs/2;  fT = linspace(0, fs/2, 800)';
plot(f1(h1), 10*log10(P1(h1)), 'Color', [0.6 0.6 0.6]); hold on;
plot(fW(h2), 10*log10(Pw1(h2)), 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
plot(fT, 10*log10(Strue(fT)), 'r--', 'LineWidth', 1.5);
xlim([0 12]); ylim([-40 12]); grid on;
xlabel('frequency (MHz)'); ylabel('PSD (dB re peak)');
legend('periodogram N = 1024', 'Welch M = 128, 50 %', 'true S(f)', 'Location', 'northeast');
title('One record: raw periodogram vs Welch');
subplot(1, 2, 2);
edges = 0:0.25:6;  cnt = histc(Rkeep, edges);  cnt = cnt(1:end-1);
bar(edges(1:end-1) + 0.125, cnt/(numel(Rkeep)*0.25), 1, 'FaceColor', [0.7 0.8 1]); hold on;
xx = linspace(0, 6, 200);  plot(xx, exp(-xx), 'r', 'LineWidth', 1.5);
xlabel('P(f) / S(f)'); ylabel('probability density'); grid on;
title('Periodogram bins are exponential (\chi^2_2/2)');
legend('N = 1024, in-band bins', 'e^{-x}');
