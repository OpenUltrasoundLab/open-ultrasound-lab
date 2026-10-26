%% m09_fir_vs_iir_pulse.m
% Module 09 - The same RF band-pass spec met by an equiripple FIR and by a
% Chebyshev I IIR: order, delay, group-delay dispersion and pulse distortion.
% Spec: pass 3-7 MHz (0.5 dB), stop <= 2 and >= 8 MHz (40 dB), fs = 40 MHz.
%  FIR: order 72 equiripple by Lawson IRLS (linear phase, delay 36 samples)
%  IIR: Chebyshev I prototype N = 8 -> band-pass order 16, prewarped bilinear, biquads
% Pulse: Gaussian-enveloped 5 MHz, 40 % (-6 dB) fractional bandwidth.
% What to change: f0, bw, the IIR prototype order, forward-backward on/off.
% Expected output: FIR group delay 36.0 samples flat; IIR group delay ~21 samples at
% 5 MHz rising to >100 near the 3 MHz edge; rho(FIR) ~ 1.000, rho(IIR causal) clearly
% lower than rho(IIR magnitude only); forward-backward IIR aligned with zero shift.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.

clear; close all; clc;
rng(1);

fs = 40;                                    % MHz
W  = @(f) 2*pi*f/fs;                        % MHz -> rad/sample
Ap = 0.5; As = 40;
dp = (10^(Ap/20) - 1)/(10^(Ap/20) + 1); ds = 10^(-As/20);

%% FIR: type I equiripple via IRLS, order 72
M = 72; L = M/2; ng = 2000;
tot = W(2) + (W(7) - W(3)) + (pi - W(8));               % grid per band, edges included exactly
g = [linspace(0, W(2), round(ng*W(2)/tot))'; linspace(W(3), W(7), round(ng*(W(7)-W(3))/tot))'; linspace(W(8), pi, round(ng*(pi-W(8))/tot))'];
d = double(g >= W(3) - 1e-12 & g <= W(7) + 1e-12); Wt = ones(size(g)); Wt(d == 0) = dp/ds;
C = cos(g * (0:L)); u = Wt;
for it = 1:80
    a = (C .* u) \ (d .* u); e = abs(Wt .* (C*a - d)); u = u .* sqrt(e + 1e-14); u = u / sqrt(mean(u.^2));
end
h = [flipud(a(2:end))/2; a(1); a(2:end)/2];

%% IIR: Chebyshev I prototype -> band-pass -> bilinear, as biquads
N = 8; ep = sqrt(10^(Ap/10) - 1); mu = asinh(1/ep)/N; th = pi*(2*(1:N) - 1)/(2*N);
pl = -sinh(mu)*sin(th) + 1j*cosh(mu)*cos(th);              % prototype poles, passband edge 1 rad/s
W1 = 2*tan(W(3)/2); W2 = 2*tan(W(7)/2); W0 = sqrt(W1*W2); B = W2 - W1;
pb = [pl*B/2 + sqrt((pl*B/2).^2 - W0^2), pl*B/2 - sqrt((pl*B/2).^2 - W0^2)];   % s -> (s^2 + W0^2)/(B s)
pz = (2 + pb) ./ (2 - pb);                                   % bilinear (T = 1)
pu = pz(imag(pz) > 1e-12);                                  % one of each conjugate pair (16 poles -> 8)
% each biquad: zeros at z = +1 and z = -1 (analog zeros at s = 0 and s = inf)
evalH = @(w) prod((exp(-2j*w(:)) - 1) * ones(1, numel(pu)) ./ (1 - 2*real(pu(:).') .* exp(-1j*w(:)) + abs(pu(:).').^2 .* exp(-2j*w(:))), 2);
wg = linspace(0, pi, 8001)'; fg = wg/(2*pi)*fs;
Hraw = evalH(wg); gIIR = 1/max(abs(Hraw(fg >= 3 & fg <= 7)));  % Chebyshev I: passband peak = 1

%% Responses and group delays
Hfir = exp(-1j*wg*(0:M)) * h;
Hiir = gIIR * Hraw;
gdf = @(H) -gradient(unwrap(angle(H)), wg(2) - wg(1));
gdFIR = gdf(Hfir); gdIIR = gdf(Hiir);
sel = fg >= 3 & fg <= 7;
fprintf('FIR order %d: passband dev %.4f (spec %.4f), stopband %.2f dB (spec -%g)\n', M, max(abs(abs(Hfir(sel)) - 1)), dp, 20*log10(max(abs(Hfir(fg <= 2 | fg >= 8)))), As);
fprintf('IIR Chebyshev I band-pass order %d: passband loss %.3f dB, stopband %.2f dB\n', 2*N, -20*log10(min(abs(Hiir(sel)))), 20*log10(max(abs(Hiir(fg <= 2 | fg >= 8)))));
[~, i5] = min(abs(fg - 5));
fprintf('Group delay: FIR %.1f samples everywhere in 3-7 MHz (%.2f to %.2f)\n', gdFIR(i5), min(gdFIR(sel)), max(gdFIR(sel)));
fprintf('             IIR %.1f samples at 5 MHz, %.1f to %.1f over 3-7 MHz\n', gdIIR(i5), min(gdIIR(sel)), max(gdIIR(sel)));
fprintf('Multiplies per output: FIR %d (folded), IIR %d (5 per biquad)\n', M/2 + 1, 5*numel(pu));

%% Pulse through both
Ns = 1024; n0 = 200; f0 = 5; bw = 0.40;
st = sqrt(2*log(2)) / (pi*bw*f0) * fs;                      % envelope sigma (samples)
n = (0:Ns-1)'; x = exp(-(n - n0).^2/(2*st^2)) .* cos(2*pi*f0*(n - n0)/fs);
yF = filter(h, 1, x);
yI = x; for i = 1:numel(pu), yI = filter([1 0 -1], [1, -2*real(pu(i)), abs(pu(i))^2], yI); end; yI = gIIR * yI;
yB = flipud(yI);                                            % forward-backward: reverse the causal output,
for i = 1:numel(pu), yB = filter([1 0 -1], [1, -2*real(pu(i)), abs(pu(i))^2], yB); end   % filter again,
yB = flipud(gIIR * yB);                                     % and reverse back: |H|^2, zero phase
% magnitude-only version (|H| applied with zero phase) via FFT
X = fft(x); kf = (0:Ns-1)'; wk = 2*pi*min(kf, Ns - kf)/Ns; Hm = abs(gIIR * evalH(wk)); yM = real(ifft(X .* Hm));
% envelope by FFT analytic signal
anal = @(v) ifft(fft(v) .* [1; 2*ones(Ns/2 - 1, 1); 1; zeros(Ns/2 - 1, 1)]);
env = @(v) abs(anal(v));
% correlation coefficient after best integer alignment (circular is safe: signals are zero-padded)
rho = @(y) max(real(ifft(fft(y) .* conj(fft(x))))) / (norm(y) * norm(x));
[~, kF] = max(env(yF)); [~, kI] = max(env(yI)); [~, kB] = max(env(yB)); [~, kx] = max(env(x));
fprintf('\nPulse f0 = %g MHz, %d %% bandwidth:\n', f0, round(100*bw));
fprintf('  envelope peak delay: FIR %d, IIR %d, IIR forward-backward %d samples\n', kF - kx, kI - kx, kB - kx);
fprintf('  rho: FIR %.3f | IIR causal %.3f | IIR magnitude only %.3f | IIR forward-backward %.3f\n', rho(yF), rho(yI), rho(yM), rho(yB));

%% Figures
figure('Name', 'FIR vs IIR band-pass');
subplot(3,1,1); plot(fg, 20*log10(abs(Hfir) + 1e-12), 'b', fg, 20*log10(abs(Hiir) + 1e-12), 'r', 'LineWidth', 1.1); hold on;
plot([0 2], [-As -As], 'k--', [8 20], [-As -As], 'k--'); ylim([-80 5]); grid on; xlabel('frequency (MHz)'); ylabel('|H| (dB)'); title('Equiripple FIR (blue) vs Chebyshev I (red)');
gdFp = gdFIR; gdFp(fg < 2.5 | fg > 7.5) = NaN;              % FIR phase jumps at its stopband zeros
subplot(3,1,2); plot(fg, gdFp, 'b', fg, gdIIR, 'r', 'LineWidth', 1.1); xlim([1 9]); ylim([0 150]); grid on; xlabel('frequency (MHz)'); ylabel('\tau_g (samples)'); title('Group delay');
subplot(3,1,3); t = (n - n0)/fs; plot(t, circshift(x, 36), 'k', t, yF, 'b', t, yI, 'r', t, yB, 'm--', 'LineWidth', 1); xlim([-1 3]); grid on;
xlabel('time (\mus)'); ylabel('amplitude'); title('Input (shifted 36 samples), FIR, IIR causal, IIR forward-backward');
