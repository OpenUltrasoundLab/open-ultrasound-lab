%% m12_frequency_crlb.m
% Purpose : frequency of a tone in noise (slow-time Doppler is this problem).
%   Part 1: complex tone y[n] = A exp(j(w n + phi)) + w[n], n = 0..N-1, in
%           circular complex Gaussian noise of variance sigma^2. The MLE is
%           the periodogram peak (zero-padded FFT, then a fine grid and a
%           parabola). Compared with the CRLB
%               var(w_hat) >= 6 / (SNR N (N^2 - 1)),  SNR = A^2/sigma^2,
%           and with the lag-one (Kasai) estimator. Shows the threshold.
%   Part 2: real sinusoid x[n] = A cos(2 pi f n + phi) + w[n] (Kay, Vol. I,
%           Example 3.14): var(f_hat) >= 12 / ((2 pi)^2 eta N (N^2 - 1)),
%           eta = A^2/(2 sigma^2). The exact MLE (projection on cos/sin)
%           reaches it; the plain periodogram peak does not at short N.
% Change  : N, the SNR grid, f0/PRF for the velocity conversion.
% Expect  : N = 32: MSE/CRLB ~1.0 for SNR >= -2 dB (N*SNR >= 13 dB); ~180x
%           at -4 dB; hundreds of times below. Kasai approaches
%           N(N+1)/(6(N-1)) = 5.7x the CRLB at high SNR and is far worse at
%           moderate SNR. Part 2: exact MLE ~1.0x; periodogram peak 1.4x,
%           4.3x, 35x at eta = 10, 20, 30 dB (bias from the -f image).
rng(2);

%% Part 1: complex tone, MLE and Kasai vs CRLB
N = 32;                       % samples (ensemble length)
snrdB = -10:2:24;             % per-sample SNR A^2/sigma^2 (dB)
nTrial = 1000;
L = 2^nextpow2(16*N);         % zero-padded FFT length
n = (0:N-1)';
nFine = 41;                   % fine grid points across +-1 FFT bin
mseMLE = zeros(size(snrdB)); mseKas = mseMLE; outFrac = mseMLE; crlb = mseMLE;
for k = 1:numel(snrdB)
    snr = 10^(snrdB(k)/10);
    w0 = (rand(1, nTrial) - 0.5) * 2*pi;          % true frequencies (rad/sample)
    ph = rand(1, nTrial) * 2*pi;
    Y = exp(1j*(n*w0 + ones(N, 1)*ph)) + sqrt(1/(2*snr)) * (randn(N, nTrial) + 1j*randn(N, nTrial));
    P = abs(fft(Y, L)).^2;
    [~, kmax] = max(P);
    wc = 2*pi*(kmax - 1)/L;                       % coarse peak
    errM = zeros(1, nTrial);
    for t = 1:nTrial
        wg = wc(t) + linspace(-2*pi/L, 2*pi/L, nFine);   % fine grid
        J = abs(exp(-1j*n*wg).' * Y(:, t)).^2;            % periodogram on the grid
        [~, j] = max(J); j = min(max(j, 2), nFine - 1);
        dlt = 0.5*(J(j-1) - J(j+1)) / (J(j-1) - 2*J(j) + J(j+1));
        wHat = wg(j) + dlt*(wg(2) - wg(1));
        errM(t) = angle(exp(1j*(wHat - w0(t))));  % wrap to (-pi, pi]
    end
    errK = angle(sum(Y(2:end, :) .* conj(Y(1:end-1, :)), 1) .* exp(-1j*w0));
    crlb(k) = 6 / (snr * N * (N^2 - 1));
    mseMLE(k) = mean(errM.^2); mseKas(k) = mean(errK.^2);
    outFrac(k) = mean(abs(errM) > 0.5*2*pi/N);    % beyond half a mainlobe
end
fprintf('Part 1: complex tone, N = %d, %d trials per SNR\n', N, nTrial);
fprintf('  SNR(dB)  N*SNR(dB)  MLE/CRLB  Kasai/CRLB  outliers\n');
for k = 1:numel(snrdB)
    fprintf('  %6.0f  %8.1f  %9.2f  %9.2f  %8.3f\n', snrdB(k), snrdB(k) + 10*log10(N), ...
        mseMLE(k)/crlb(k), mseKas(k)/crlb(k), outFrac(k));
end
% Doppler conversion at one operating point
c = 1540; f0 = 5e6; PRF = 5e3; M = 16; snr = 10;   % 10 dB per pulse
sw = sqrt(6/(snr*M*(M^2 - 1)));
fprintf('  Doppler example: M = %d, SNR = 10 dB: sigma_w = %.4f rad/sample, sigma_fD = %.2f Hz, sigma_v = %.2f mm/s\n', ...
    M, sw, sw*PRF/(2*pi), c*sw*PRF/(2*pi)/(2*f0)*1e3);

%% Part 2: real sinusoid, exact MLE vs periodogram peak vs Kay's CRLB
N2 = 32; n2 = (0:N2-1)'; etadB = [10 20 30]; nT2 = 2000; L2 = 512;
fprintf('\nPart 2: real sinusoid, N = %d, f0 ~ U(0.1, 0.4) cycles/sample\n', N2);
for k = 1:numel(etadB)
    eta = 10^(etadB(k)/10); s2 = 1/(2*eta);              % A = 1
    f0t = 0.1 + 0.3*rand(1, nT2); ph = 2*pi*rand(1, nT2);
    X = cos(2*pi*n2*f0t + ones(N2, 1)*ph) + sqrt(s2)*randn(N2, nT2);
    eP = zeros(1, nT2); eX = eP;
    for t = 1:nT2
        x = X(:, t);
        Pq = abs(fft(x, L2)).^2; [~, kk] = max(Pq(2:L2/2)); kk = kk + 1;
        fg = (kk - 1)/L2 + linspace(-1/L2, 1/L2, nFine);
        E = exp(-1j*2*pi*n2*fg);
        Jp = abs(E.' * x).^2;                           % periodogram on the grid
        C = real(E); S = -imag(E);                       % cos and sin columns
        xc = C.'*x; xs = S.'*x; cc = sum(C.^2).'; ss = sum(S.^2).'; cs = sum(C.*S).';
        Jx = (ss.*xc.^2 - 2*cs.*xc.*xs + cc.*xs.^2) ./ (cc.*ss - cs.^2);   % ||P_H x||^2
        [~, j] = max(Jp); j = min(max(j, 2), nFine - 1);
        eP(t) = fg(j) + 0.5*(Jp(j-1) - Jp(j+1))/(Jp(j-1) - 2*Jp(j) + Jp(j+1))*(fg(2) - fg(1)) - f0t(t);
        [~, j] = max(Jx); j = min(max(j, 2), nFine - 1);
        eX(t) = fg(j) + 0.5*(Jx(j-1) - Jx(j+1))/(Jx(j-1) - 2*Jx(j) + Jx(j+1))*(fg(2) - fg(1)) - f0t(t);
    end
    crlbR = 12 / ((2*pi)^2 * eta * N2 * (N2^2 - 1));
    fprintf('  eta = %2d dB: exact MLE var/CRLB = %.2f, periodogram peak MSE/CRLB = %.2f\n', ...
        etadB(k), var(eX)/crlbR, mean(eP.^2)/crlbR);
end

%% Figure
figure('Position', [80 80 900 380]);
subplot(1, 2, 1);
semilogy(snrdB, sqrt(mseMLE), 'bo-', snrdB, sqrt(mseKas), 's-', 'Color', [0.8 0.5 0]); hold on;
semilogy(snrdB, sqrt(crlb), 'k-', 'LineWidth', 1.4);
semilogy(snrdB, pi/sqrt(3)*ones(size(snrdB)), 'k:');
xlabel('per-sample SNR (dB)'); ylabel('RMS error (rad/sample)');
title(sprintf('Complex tone, N = %d', N)); legend('periodogram MLE', 'Kasai (lag one)', 'CRLB', 'uniform guess pi/sqrt(3)');
subplot(1, 2, 2);
plot(snrdB + 10*log10(N), outFrac, 'bo-'); grid on;
xlabel('integrated SNR N\cdotSNR (dB)'); ylabel('outlier fraction');
title('Threshold: outliers beyond half a mainlobe');
