%% m07_speckle_autocorr_psd.m  --  Wiener-Khinchin, speckle correlation length and N_eff
% Purpose : simulate RF speckle lines as white scatterer reflectivity filtered by a
%           Gaussian pulse (f0 = 5 MHz, 60 % -6 dB bandwidth, fs = 50 MHz). Check that
%           (1) the periodogram is exactly the DFT of the biased autocorrelation
%           (Wiener-Khinchin for finite records), (2) the RF autocorrelation envelope has
%           FWHM 1.25/B and the intensity autocovariance FWHM 0.88/B, (3) the periodogram
%           does not get less noisy with a longer record, and (4) the effective number of
%           independent samples in a 5 mm window depends on the statistic: the speckle-cell
%           count T*B/1.25, the intensity mean T*B/0.94, the RF correlation coefficient
%           T*B/0.47.
% Change  : f0, frac_bw, win_mm (analysis window), rho (for the correlation estimate).
% Expect  : max |periodogram - DFT(R_biased)| ~ 1e-13 (relative); tau_FWHM*B = 1.25 (RF),
%           0.88 (intensity); axial cell 0.32 mm; periodogram std/mean ~ 1.0 for 256 and
%           4096 samples; N_eff for 5 mm = 15.6 (cells), ~21 (intensity mean, theory 20.8),
%           ~40 (rho_hat at rho = 0.9, theory 41.5), out of 325 samples.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);

%% 1. Pulse spectrum and speckle lines
c0 = 1540;  fs = 50e6;  f0 = 5e6;  frac_bw = 0.6;     % m/s, Hz, Hz, -6 dB fractional bandwidth
B   = frac_bw * f0;                                   % -6 dB (amplitude) bandwidth (Hz)
sgt = sqrt(2*log(2)) / (pi*B);                        % Gaussian envelope std of the pulse (s)
Ns  = 8192;  Nl = 200;                                % samples per line, number of lines
f   = [0:Ns/2-1, -Ns/2:-1] * fs / Ns;                 % FFT frequency grid (Hz)
Pf  = exp(-2*pi^2*sgt^2*(f - f0).^2) + exp(-2*pi^2*sgt^2*(f + f0).^2);  % pulse spectrum (amplitude)
rfl = @(n) real(ifft(fft(randn(Ns, n)) .* Pf.'));     % n lines of white reflectivity * pulse (circular)
X   = rfl(Nl);  X = X / sqrt(mean(X(:).^2));          % Ns x Nl RF, unit power
fprintf('pulse sigma_t %.1f ns, B = %.2f MHz, 1/B = %.0f ns\n', 1e9*sgt, B/1e6, 1e9/B);

%% 2. Wiener-Khinchin for a finite record: periodogram = DFT of the biased autocorrelation
x   = X(1:1024, 1);  Nx = numel(x);
Rb  = zeros(2*Nx-1, 1);                               % biased R(k), k = -(Nx-1)..(Nx-1)
for k = 0:Nx-1
    Rb(Nx+k) = sum(x(1+k:Nx) .* x(1:Nx-k)) / Nx;  Rb(Nx-k) = Rb(Nx+k);
end
Lf  = 2*Nx - 1;  kk = (0:Lf-1).';
Sper = abs(fft(x, Lf)).^2 / Nx;                        % periodogram on a 2Nx-1 point grid
Swk  = real(exp(1i*2*pi*kk*(Nx-1)/Lf) .* fft(Rb));     % DFT of R(k), index shifted to k = -(Nx-1)
fprintf('max |periodogram - DFT(R_biased)| / max = %.1e\n', max(abs(Sper - Swk)) / max(Sper));

%% 3. Correlation lengths: RF envelope (1.25/B) and intensity (0.88/B)
Xf  = fft(X);  Xf(Ns/2+2:end, :) = 0;  Xf(2:Ns/2, :) = 2*Xf(2:Ns/2, :);
Z   = ifft(Xf);                                       % analytic signal (IQ at RF), Ns x Nl
Ra  = mean(ifft(abs(fft(Z, 2*Ns)).^2), 2);            % averaged complex autocorrelation (biased)
env = abs(Ra) / abs(Ra(1));                           % RF autocorrelation envelope vs lag
Ia  = abs(Z).^2;  Ia = Ia - mean(Ia(:));              % intensity fluctuation
RI  = mean(ifft(abs(fft(Ia, 2*Ns)).^2), 2);  RI = real(RI) / real(RI(1));
fwhm = @(r) 2 * interp1(r(find(r < 0.5, 1) + (-1:0)), (find(r < 0.5, 1) + (-1:0)) - 1, 0.5) / fs;  % s
tR = fwhm(env);  tI = fwhm(RI);
fprintf('RF  autocorrelation FWHM %.0f ns, x B = %.3f (theory 4*sqrt(2)*ln2/pi = %.3f)\n', 1e9*tR, tR*B, 4*sqrt(2)*log(2)/pi);
fprintf('intensity autocov.  FWHM %.0f ns, x B = %.3f (theory /sqrt(2) = %.3f)\n', 1e9*tI, tI*B, 4*log(2)/pi);
fprintf('axial speckle cell c*tau/2 = %.3f mm (RF), %.3f mm (intensity)\n', 1e3*c0*tR/2, 1e3*c0*tI/2);

%% 4. The periodogram does not converge: each bin is ~exponential, std/mean ~ 1
for Nw = [256 4096]
    P = abs(fft(X(1:Nw, :))).^2 / Nw;                  % Nw x Nl periodograms
    b = round([4.5 5 5.5] * 1e6 / fs * Nw) + 1;
    fprintf('periodogram, %4d samples: std/mean across lines at 4.5/5/5.5 MHz = %.2f %.2f %.2f\n', ...
        Nw, std(P(b, :), 0, 2) ./ mean(P(b, :), 2));
end

%% 5. Effective number of independent samples in a 5 mm window
win_mm = 5;  W = round(2*win_mm*1e-3/c0*fs);  T = W/fs;  % window: samples, duration (s)
nw = floor(Ns / W);
Iw = reshape(abs(Z(1:nw*W, :)).^2, W, []);            % intensity, one column per window
Neff_I = var(abs(Z(:)).^2) / var(mean(Iw, 1));          % sigma_I^2 / var(window mean)
rho = 0.9;  Nn = rfl(Nl);  Nn = Nn / sqrt(mean(Nn(:).^2));   % independent speckle, unit power
Y = rho*X + sqrt(1 - rho^2) * Nn;                     % second signal, correlation rho with X
Xw = reshape(X(1:nw*W, :), W, []);  Yw = reshape(Y(1:nw*W, :), W, []);
Xw = Xw - mean(Xw, 1);  Yw = Yw - mean(Yw, 1);
r  = sum(Xw .* Yw, 1) ./ sqrt(sum(Xw.^2, 1) .* sum(Yw.^2, 1));
Neff_r = ((1 - rho^2) / std(r))^2;                    % from std(rho_hat) = (1 - rho^2)/sqrt(N_eff)
fprintf('5 mm window = %d samples, T*B = %.2f\n', W, T*B);
fprintf('N_eff: speckle cells T/tau_FWHM %.1f | intensity mean %.1f (theory T/(sqrt(2pi) sigma_t) = %.1f) | rho_hat %.1f (theory 2T/(sqrt(2pi) sigma_t) = %.1f)\n', ...
    T*B/1.248, Neff_I, T/(sqrt(2*pi)*sgt), Neff_r, 2*T/(sqrt(2*pi)*sgt));
fprintf('std(rho_hat) at rho = %.1f: measured %.4f | from cells %.4f | from samples %.4f\n', ...
    rho, std(r), (1 - rho^2)/sqrt(T*B/1.248 - 1), (1 - rho^2)/sqrt(W - 1));

%% 6. Figures
figure('Position', [100 100 1100 380]);
subplot(1,3,1);
t = (0:399) / fs * 1e6;  plot(t, X(1:400, 1), 'b-'); grid on;
xlabel('time (\mus)'); ylabel('RF (\sigma units)'); title('speckle RF line (first 8 \mus)');
subplot(1,3,2);
lagk = -80:80;  lt = lagk / fs * 1e6;
Rr  = real(Ra);  Rr = Rr / Rr(1);
idx = mod(lagk, 2*Ns) + 1;
plot(lt, Rr(idx), 'b-', 'LineWidth', 1.4); hold on;
plot(lt, exp(-(lagk/fs).^2 / (4*sgt^2)) .* cos(2*pi*f0*lagk/fs), 'r--', 'LineWidth', 1.2);
plot(lt, env(idx), 'k-');  plot(lt, RI(idx), '-', 'Color', [0 0.6 0.4], 'LineWidth', 1.4);
grid on; xlim([-1.6 1.6]); ylim([-1.05 1.1]);
xlabel('lag \tau (\mus)'); ylabel('normalized correlation');
title('RF autocorr. (blue), theory (red), envelope (black), intensity (green)');
subplot(1,3,3);
Nw = 256;
Xs = reshape(X(1:Nw*floor(Ns/Nw), :), Nw, []);         % non-overlapping 256-sample segments
hw = 0.5 - 0.5*cos(2*pi*(0:Nw-1).'/Nw);                 % Hann window (written out)
Pw = mean(abs(fft(Xs .* hw)).^2, 2);  fw = (0:Nw/2) * fs / Nw;
Pth = Pf(1:Ns/2+1).^2;  fth = f(1:Ns/2+1);
plot(fw/1e6, 10*log10(Pw(1:Nw/2+1) / max(Pw)), 'b-', 'LineWidth', 1.4); hold on;
plot(fth/1e6, 10*log10(Pth / max(Pth)), 'r--', 'LineWidth', 1.2);
grid on; xlim([0 15]); ylim([-40 3]);
xlabel('frequency (MHz)'); ylabel('PSD (dB re peak)'); title('averaged periodogram vs |P(f)|^2');
