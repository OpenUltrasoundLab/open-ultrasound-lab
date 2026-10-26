%% m27_ar_burg.m  -  parametric (AR) spectra for short gates: Burg vs Yule-Walker vs periodogram
% Purpose : (1) two tones 0.4 MHz apart in a 64-sample gate (1.6 us at 40 MHz,
%               1.23 mm of depth): the periodogram cannot separate them
%               (Rayleigh resolution fs/N = 0.625 MHz); Burg AR(p) can once p
%               is large enough; Yule-Walker (biased autocorrelation = implicit
%               triangular window) needs a much higher order;
%           (2) AIC order selection: N ln(E_p) + 2p;
%           (3) on a short speckle gate, the spectral centroid from Burg AR vs a
%               Hann periodogram: both unbiased, AR about 20 % lower std, and
%               the AR order hardly matters because the one-sided centroid is
%               fs/4 - (2 fs/pi^2) sum_{odd l} rho(l)/l^2, dominated by lag 1,
%               and every Burg model of order >= l shares rho(0..l).
% Change  : N, tone separation, SNR, model order p.
% Expect  : P(resolved) Burg 0.37 / 0.66 / 0.75 at p = 12 / 16 / 24;
%           Yule-Walker 0.00 / 0.09 / 0.29; Hann periodogram 0.02; AIC picks
%           p = 16 (median). Speckle centroid std 0.36 MHz (periodogram) vs
%           0.29 MHz (Burg, any order).
clear; close all; rng(1);  clear x1;

%% Parameters
fs = 40;  N = 64;  f1 = 4.8;  f2 = 5.2;     % MHz, samples
snr_dB = 20;                                % per tone: (A^2/2)/sigma^2 with A = 1
ntr = 200;  orders = [4 8 12 16 24];
fg = linspace(3, 7, 801)';                  % evaluation grid (MHz)
E1 = exp(-1i*2*pi*fg/fs*(0:max(orders)));   % e^{-j 2 pi f k / fs}
arpsd = @(a, E) E./abs(E1(:, 1:numel(a))*a).^2;
% Burg recursion is written out inline below (scripts cannot hold local functions portably)

%% 1. Resolution test
fprintf('Two tones %.1f and %.1f MHz, N = %d (%.2f us), SNR %d dB per tone, %d trials\n', f1, f2, N, N/fs, snr_dB, ntr);
fprintf('Periodogram Rayleigh resolution fs/N = %.3f MHz > separation %.1f MHz\n', fs/N, f2 - f1);
resB = zeros(numel(orders), 2);  resY = resB;  resP = 0;  aicSel = zeros(ntr, 1);
n = (0:N-1)';  sig = sqrt(0.5*10^(-snr_dB/10));
for t = 1:ntr
  x = cos(2*pi*f1/fs*n + 2*pi*rand) + cos(2*pi*f2/fs*n + 2*pi*rand) + sig*randn(N, 1);
  % --- Burg up to max order, storing every intermediate model
  pmax = max(orders);  ef = x;  eb = x;  a = 1;  Ep = zeros(pmax+1, 1);  Ep(1) = x'*x/N;  Aall = cell(pmax, 1);
  for m = 0:pmax-1
    ff = ef(m+2:N);  bb = eb(m+1:N-1);
    k = -2*(bb'*ff)/(ff'*ff + bb'*bb);
    a = [a; 0] + k*flipud([a; 0]);
    ef(m+2:N) = ff + k*bb;  eb(m+2:N) = bb + k*ff;
    Ep(m+2) = Ep(m+1)*(1 - k^2);  Aall{m+1} = a;
  end
  aic = N*log(Ep) + 2*(0:pmax)';  [~, im] = min(aic);  aicSel(t) = im - 1;
  % --- Yule-Walker with the biased autocorrelation estimate
  r = zeros(pmax+1, 1);  for l = 0:pmax, r(l+1) = x(1:N-l)'*x(1+l:N)/N; end
  for io = 1:numel(orders)
    p = orders(io);
    for meth = 1:2
      if meth == 1, ap = Aall{p}; Pp = arpsd(ap, Ep(p+1));
      else, ay = [1; -toeplitz(r(1:p))\r(2:p+1)]; Pp = arpsd(ay, r(1) + r(2:p+1)'*ay(2:end)); end
      pk = find(Pp(2:end-1) > Pp(1:end-2) & Pp(2:end-1) >= Pp(3:end) & Pp(2:end-1) > 0.01*max(Pp)) + 1;
      ok = numel(pk) == 2 && all(abs(sort(fg(pk)) - [f1; f2]) < 0.15);
      if meth == 1 && p == 16 && ok && ~exist('x1', 'var'), x1 = x; A1 = Aall; E1p = Ep; end
      if meth == 1, resB(io, :) = resB(io, :) + [ok, numel(pk) > 2]; else, resY(io, :) = resY(io, :) + [ok, numel(pk) > 2]; end
    end
  end
  Px = abs(exp(-1i*2*pi*fg/fs*n')*(x.*(0.5 - 0.5*cos(2*pi*(n + 0.5)/N)))).^2;
  pk = find(Px(2:end-1) > Px(1:end-2) & Px(2:end-1) >= Px(3:end) & Px(2:end-1) > 0.01*max(Px)) + 1;
  resP = resP + (numel(pk) == 2 && all(abs(sort(fg(pk)) - [f1; f2]) < 0.15));
end
fprintf('\n order p   Burg: P(resolved)  P(>2 peaks)   Yule-Walker: P(resolved)\n');
for io = 1:numel(orders)
  fprintf('   %3d        %5.3f          %5.3f              %5.3f\n', orders(io), resB(io, 1)/ntr, resB(io, 2)/ntr, resY(io, 1)/ntr);
end
fprintf(' Hann periodogram: P(resolved) = %.3f\n', resP/ntr);
qs = sort(aicSel);
fprintf(' AIC-selected Burg order: median %d, 10-90 %% range %d-%d\n', median(aicSel), qs(round(0.1*ntr)), qs(round(0.9*ntr)));

%% 3. Short speckle gate: centroid from Burg AR(p) vs Hann periodogram
f0 = 5;  sigf = 0.6*f0/(2*sqrt(2*log(2)));  Nf = 1024;
fF = (0:Nf-1)'*fs/Nf;  fF(fF >= fs/2) = fF(fF >= fs/2) - fs;
H = exp(-(abs(fF) - f0).^2/(2*sigf^2));
fc = linspace(0.5, 12, 1151)';  Ec = exp(-1i*2*pi*fc/fs*(0:16));
cen = @(P) sum(fc.*P)/sum(P);
ntr2 = 500;  Cp = zeros(ntr2, 1);  Cb = zeros(ntr2, 3);  pl = [4 8 16];
wN = 0.5 - 0.5*cos(2*pi*(n + 0.5)/N);  En = exp(-1i*2*pi*fc/fs*n');
for t = 1:ntr2
  x = real(ifft(fft(randn(Nf, 1)).*H));  x = x(1:N);
  Cp(t) = cen(abs(En*(x.*wN)).^2);
  ef = x;  eb = x;  a = 1;  Ep = x'*x/N;
  for m = 0:max(pl)-1
    ff = ef(m+2:N);  bb = eb(m+1:N-1);  k = -2*(bb'*ff)/(ff'*ff + bb'*bb);
    a = [a; 0] + k*flipud([a; 0]);  ef(m+2:N) = ff + k*bb;  eb(m+2:N) = bb + k*ff;  Ep = Ep*(1 - k^2);
    j = find(pl == m + 1);  if ~isempty(j), Cb(t, j) = cen(Ep./abs(Ec(:, 1:numel(a))*a).^2); end
  end
end
fprintf('\nSpeckle gate N = %d (%.2f mm), true power-spectrum centroid %.2f MHz, %d gates:\n', N, N/fs*1.54/2, f0, ntr2);
fprintf('  Hann periodogram : mean %.3f MHz, std %.3f MHz\n', mean(Cp), std(Cp));
for j = 1:numel(pl), fprintf('  Burg AR(%2d)      : mean %.4f MHz, std %.4f MHz\n', pl(j), mean(Cb(:, j)), std(Cb(:, j))); end

%% Figures
figure('Position', [100 100 950 380]);
subplot(1, 2, 1);
Px1 = abs(exp(-1i*2*pi*fg/fs*n')*(x1.*(0.5 - 0.5*cos(2*pi*(n + 0.5)/N)))).^2;
plot(fg, 10*log10(Px1/max(Px1)), 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2); hold on;
for p = [8 16]
  Pp = arpsd(A1{p}, E1p(p+1));  plot(fg, 10*log10(Pp/max(Pp)), 'LineWidth', 1.3);
end
plot([f1 f1], [-40 3], 'k:', [f2 f2], [-40 3], 'k:');
ylim([-40 3]); grid on; xlabel('frequency (MHz)'); ylabel('dB re peak');
legend('Hann periodogram', 'Burg p = 8', 'Burg p = 16', 'true tones', 'Location', 'south');
title(sprintf('A gate Burg(16) resolves: N = %d, tones %.1f MHz apart', N, f2 - f1));
subplot(1, 2, 2);
plot(orders, resB(:, 1)/ntr, 'o-', orders, resY(:, 1)/ntr, 's-', 'LineWidth', 1.3); grid on;
xlabel('AR order p'); ylabel('P(two peaks within 0.15 MHz)'); ylim([0 1]);
legend('Burg', 'Yule-Walker (biased r)', 'Location', 'northwest');
title('Resolution probability vs order');
