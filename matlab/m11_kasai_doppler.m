%% m11_kasai_doppler.m
% Purpose : lag-one autocorrelation (Kasai) velocity estimator on a slow-time
%           IQ ensemble: (1) estimate vs true velocity with aliasing beyond
%           v_Nyq = c*PRF/(4 f0); (2) velocity std vs ensemble length M,
%           compared with a high-SNR formula and with the CRLB for a tone.
% Change  : PRF, f0, SNR, M. Doubling PRF doubles v_Nyq and, at fixed phase
%           noise, doubles the velocity std.
% Expect  : f0 = 5 MHz, PRF = 5 kHz: v_Nyq = 0.385 m/s, lambda/4 = 77 um per
%           pulse. Part 2 (one gate, constant amplitude, SNR 20 dB): std
%           follows var(dphi) = 1/(SNR (M-1)^2) + 1/(2 SNR^2 (M-1)), so it falls
%           as 1/(M-1) when signal x noise terms dominate; it equals the CRLB
%           6/(SNR M (M^2-1)) at M = 2 and 3 and is about 1.8x the CRLB in
%           variance at M = 8: at high SNR the phase sum telescopes, so the
%           estimate effectively uses only the first and last pulses.
rng(1);

%% Parameters
c   = 1540;     % speed of sound (m/s)
f0  = 5e6;      % transmit centre frequency (Hz)
PRF = 5e3;      % pulse repetition frequency (Hz)
G   = 8;        % axial gates averaged in R(1)
M   = 8;        % ensemble length (pulses)
snrdB = 20;     % per-sample SNR (signal power / noise power)
vNyq = c * PRF / (4 * f0);
fprintf('v_Nyq = c*PRF/(4 f0) = %.3f m/s; lambda/4 = %.1f um; ensemble duration = %.2f ms\n', ...
    vNyq, c/f0/4*1e6, M/PRF*1e3);

% Kasai estimate from an ensemble z (M x G): R(1) summed over pulses and gates
kasai = @(z) c * PRF / (4*pi*f0) * angle(sum(sum(z(2:end, :) .* conj(z(1:end-1, :)))));
% Simulate a G-gate ensemble: random complex speckle amplitude per gate, constant over slow time
ens = @(v, Mm, Gg, snr) bsxfun(@times, (randn(1,Gg) + 1j*randn(1,Gg))/sqrt(2), ...
      exp(1j * (0:Mm-1)' * 4*pi*f0*v/(c*PRF))) + (randn(Mm,Gg) + 1j*randn(Mm,Gg)) / sqrt(2*snr);

%% Part 1: estimate vs true velocity (aliasing)
snr = 10^(snrdB/10);
vTrue = linspace(-1.5, 1.5, 241);
vHat = zeros(size(vTrue));
for i = 1:numel(vTrue), vHat(i) = kasai(ens(vTrue(i), M, G, snr)); end
vWrap = vTrue - 2*vNyq * round(vTrue / (2*vNyq));
fprintf('Part 1: max |v_hat - wrapped v| inside the band = %.3f m/s\n', max(abs(vHat(abs(vTrue) < 0.9*vNyq) - vWrap(abs(vTrue) < 0.9*vNyq))));
v1 = 0.5;  fprintf('v = %.2f m/s (1.3 v_Nyq) reads as %.3f m/s (expected wrap %.3f)\n', v1, kasai(ens(v1, M, G, snr)), v1 - 2*vNyq);

%% Part 2: std vs ensemble length, one gate, constant amplitude (a tone in noise)
Ms = [2 3 4 6 8 12 16 24 32 48 64];  nT = 4000;  v0 = 0.1;
sdMC = zeros(size(Ms));  sdTh = sdMC;  sdCR = sdMC;
k = c * PRF / (4*pi*f0);                         % m/s per radian
for im = 1:numel(Ms)
    Mm = Ms(im);  m = (0:Mm-1)';
    ph = 2*pi*rand(1, nT);                       % random start phase per trial
    z  = exp(1j * bsxfun(@plus, m * 4*pi*f0*v0/(c*PRF), ph)) + (randn(Mm,nT) + 1j*randn(Mm,nT)) / sqrt(2*snr);
    est = k * angle(sum(z(2:end,:) .* conj(z(1:end-1,:)), 1));
    sdMC(im) = std(est);
    sdTh(im) = k * sqrt(1/(snr*(Mm-1)^2) + 1/(2*snr^2*(Mm-1)));
    sdCR(im) = k * sqrt(6 / (snr * Mm * (Mm^2 - 1)));
end
fprintf('Part 2 (SNR %d dB, one gate, %d trials):\n   M   T_ens(ms)  sd_MC(mm/s)  sd_formula  sd_CRLB  var ratio MC/CRLB\n', snrdB, nT);
for im = 1:numel(Ms)
    fprintf('  %2d   %6.2f     %8.3f    %8.3f   %7.3f    %5.2f\n', Ms(im), Ms(im)/PRF*1e3, 1e3*sdMC(im), 1e3*sdTh(im), 1e3*sdCR(im), (sdMC(im)/sdCR(im))^2);
end

%% Figure
figure;
subplot(1, 2, 1);
vW = vWrap;  vW([false, abs(diff(vWrap)) > vNyq]) = NaN;   % break the sawtooth at the wraps
plot(vTrue, vHat, '.', vTrue, vW, '-'); hold on;
plot([vNyq vNyq], [-1.5 1.5], 'k:', [-vNyq -vNyq], [-1.5 1.5], 'k:'); hold off;
xlabel('true axial velocity (m/s)'); ylabel('Kasai estimate (m/s)'); grid on;
title(sprintf('Aliasing at v_{Nyq} = %.3f m/s', vNyq)); axis([-1.5 1.5 -0.5 0.5]);
subplot(1, 2, 2);
loglog(Ms, 1e3*sdMC, 'o', Ms, 1e3*sdTh, '-', Ms, 1e3*sdCR, 'k--');
xlabel('ensemble length M (pulses)'); ylabel('std of v (mm/s)'); grid on;
legend('Monte Carlo', 'high-SNR formula', 'CRLB (tone)', 'Location', 'southwest');
title(sprintf('Std vs M, SNR = %d dB', snrdB));
