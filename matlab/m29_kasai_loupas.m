%% m29_kasai_loupas.m
% Purpose : Kasai (1D, lag-one slow-time) versus Loupas (2D) autocorrelation
%           velocity estimators on simulated IQ ensembles of moving speckle.
%           The received pulse centre frequency fc is lower than the
%           demodulation frequency fdem (attenuation downshift). Kasai uses
%           fdem in the Doppler equation and is biased by fc/fdem; Loupas
%           estimates the mean RF frequency from the lag-one correlation along
%           depth and removes the bias. Also: the variance (spectral width)
%           estimate sigma^2 = (2/T^2)(1 - |R(1)|/R(0)) and its noise bias.
% Change  : fc, N (ensemble), M (range samples), snrdB, sigv (velocity spread).
% Expect  : fdem = 5 MHz, fc = 4.5 MHz, PRF = 4 kHz: v_Nyq = 0.308 m/s.
%           For v = 0.15 m/s Kasai reads about 0.9 x v (-10 %), Loupas about
%           1.00 x v. Std falls as N grows (roughly 1/sqrt(N) with velocity
%           spread). Noise inflates the variance estimate: with SNR s the
%           ratio |R(1)|/R(0) is multiplied by s/(1+s).
rng(1);

%% Parameters
c     = 1540;      % speed of sound (m/s)
fdem  = 5e6;       % demodulation (nominal transmit) frequency (Hz)
fc    = 4.5e6;     % received centre frequency after attenuation (Hz)
B     = 0.6;       % -6 dB fractional bandwidth of the received pulse
fs    = 20e6;      % IQ sampling rate along depth (Hz)
PRF   = 4e3;       % pulse repetition frequency (Hz)
T     = 1/PRF;
vNyq  = c*PRF/(4*fdem);
fprintf('v_Nyq (with fdem) = %.3f m/s, lambda = %.3f mm\n', vNyq, c/fdem*1e3);

% Gaussian envelope with -6 dB bandwidth B*fc: |G(f)| = exp(-(f/sf)^2/2), sf from -6 dB half width
sf  = B*fc/2 / sqrt(2*log(2));          % spectral std (Hz) of amplitude spectrum
st  = 1/(2*pi*sf);                      % envelope std in time (s)
% IQ sample of scatterers at delays tau (K x 1) with amplitudes a, observed at fast time t (1 x M):
iq = @(t, tau, a) a.' * (exp(-((t - tau).^2)/(2*st^2)) .* exp(2j*pi*(fc - fdem)*(t - tau)) .* exp(-2j*pi*fdem*tau));
% (note: t - tau uses implicit expansion, K x M)

kasai  = @(Z) c*PRF/(4*pi*fdem) * angle(sum(sum(Z(:, 2:end) .* conj(Z(:, 1:end-1)))));
loupas = @(Z) c/2 * (PRF/(2*pi) * angle(sum(sum(Z(:, 2:end) .* conj(Z(:, 1:end-1)))))) / ...
              (fdem + fs/(2*pi) * angle(sum(sum(Z(2:end, :) .* conj(Z(1:end-1, :))))));

%% Simulation loop
z0   = 30e-3;  t0 = 2*z0/c;                 % gate centre depth and its delay
dens = 20e3;                                % scatterers per metre (20 per mm)
nT   = 400;                                 % trials per setting
cases = [4 8 16];  v = 0.15;  sigv = 0.02;  M = 8;  snrdB = 20;
fprintf('v = %.3f m/s, sigma_v = %.3f m/s, M = %d range samples, SNR = %d dB, %d trials\n', v, sigv, M, snrdB, nT);
fprintf('   N   Kasai mean  bias    std     Loupas mean  bias    std     (m/s, %%)\n');
res = zeros(numel(cases), 4);
for ic = 1:numel(cases)
    N = cases(ic);  vk = zeros(nT,1);  vl = vk;
    t = t0 + ((0:M-1) - (M-1)/2)/fs;                 % fast-time samples (1 x M)
    Lw = (M/fs + 8*st) * c/2 + N*T*(abs(v) + 3*sigv);  % spatial window (m)
    for it = 1:nT
        K = max(1, round(dens*Lw));
        x0 = z0 + (rand(K,1) - 0.5)*Lw;              % initial depths (m)
        vk_s = v + sigv*randn(K,1);                  % scatterer velocities toward the probe
        a  = (randn(K,1) + 1j*randn(K,1))/sqrt(2);
        Z  = zeros(M, N);
        for n = 1:N
            tau = 2*(x0 - vk_s*(n-1)*T)/c;           % moving toward the probe: delay shrinks
            Z(:, n) = iq(t, tau, a).';
        end
        ps = mean(abs(Z(:)).^2);
        Z  = Z + sqrt(ps/10^(snrdB/10)) * (randn(M,N) + 1j*randn(M,N))/sqrt(2);
        vk(it) = kasai(Z);  vl(it) = loupas(Z);
    end
    res(ic, :) = [mean(vk) std(vk) mean(vl) std(vl)];
    fprintf('  %2d    %.4f   %+5.1f%%  %.4f    %.4f    %+5.1f%%  %.4f\n', N, mean(vk), 100*(mean(vk)/v - 1), std(vk), ...
        mean(vl), 100*(mean(vl)/v - 1), std(vl));
end
fprintf('Expected Kasai scale fc/fdem = %.3f (bias %+.1f %%)\n', fc/fdem, 100*(fc/fdem - 1));

%% Variance estimator: pure slow-time Gaussian-spectrum signal, spectral std sigf, with noise
sigf = 150;  Nv = 16;  nV = 4000;  m = (0:Nv-1)';
fprintf('Variance estimate (Gaussian spectrum, sigma_f = %d Hz, N = %d):\n  SNR(dB)   sigma_f est (Hz)   predicted with noise bias\n', sigf, Nv);
% Correlated Gaussian process by spectral shaping of white noise (long FFT, then crop)
Lf = 4096;  ff = ((0:Lf-1) - Lf/2)'/Lf*PRF;  H = exp(-ff.^2/(4*sigf^2));  H = ifftshift(H);
for snrdB2 = [40 20 10 0]
    s2 = zeros(nV,1);  snr2 = 10^(snrdB2/10);
    for it = 1:nV/40
        W = (randn(Lf,40) + 1j*randn(Lf,40))/sqrt(2);
        X = ifft(bsxfun(@times, fft(W), H));  X = X(1:Nv, :);  X = X / sqrt(mean(abs(X(:)).^2));
        X = X + (randn(Nv,40) + 1j*randn(Nv,40))/sqrt(2*snr2);
        for j = 1:40
            z = X(:, j);  R1 = sum(z(2:end).*conj(z(1:end-1)));
            R0 = (sum(abs(z(2:end)).^2) + sum(abs(z(1:end-1)).^2))/2;   % symmetric lag-0 (less finite-N bias)
            s2((it-1)*40 + j) = 2/T^2 * (1 - abs(R1)/R0);
        end
    end
    rho = snr2/(1 + snr2) * exp(-(2*pi*sigf*T)^2/2);
    fprintf('   %3d        %6.1f              %6.1f\n', snrdB2, sqrt(max(mean(s2),0))/(2*pi), sqrt(2/T^2*(1 - rho))/(2*pi));
end

%% Aliasing: Kasai reading across velocities
vv = linspace(-2.5*vNyq, 2.5*vNyq, 121);  vr = zeros(size(vv));  N = 8;
t = t0 + ((0:M-1) - (M-1)/2)/fs;
for i = 1:numel(vv)
    K = round(dens*(M/fs*c/2 + 8*st*c/2 + 0.002));  x0 = z0 + (rand(K,1)-0.5)*(M/fs*c/2 + 8*st*c/2 + 0.002);
    a = (randn(K,1) + 1j*randn(K,1))/sqrt(2);  Z = zeros(M, N);
    for n = 1:N, Z(:, n) = iq(t, 2*(x0 - vv(i)*(n-1)*T)/c, a).'; end
    vr(i) = kasai(Z);
end

%% Figures
figure('Position', [50 50 1000 360]);
subplot(1,2,1); errorbar(cases - 0.15, res(:,1), res(:,2), 'o-'); hold on;
errorbar(cases + 0.15, res(:,3), res(:,4), 's-'); plot([2 18], [v v], 'k--');
xlabel('ensemble length N'); ylabel('velocity estimate (m/s)'); legend('Kasai (f_{dem})', 'Loupas 2D', 'true v', 'Location', 'southeast');
title(sprintf('f_c = %.1f MHz, f_{dem} = %.1f MHz', fc/1e6, fdem/1e6)); grid on; xlim([2 18]);
subplot(1,2,2); plot(vv, vr, '.', vv, vv*fc/fdem, 'k:'); xlabel('true velocity (m/s)'); ylabel('Kasai estimate (m/s)');
title('Aliasing: estimates wrap at \pm v_{Nyq}'); grid on;
