%% m11_jitter_vs_crlb.m
% Purpose : Monte Carlo jitter of two delay estimators versus SNR, compared
%           with the Cramer-Rao bound (Carter coherence integral for this
%           spectrum) and the Walker-Trahey closed form; outlier fraction
%           shows the threshold effect (jumps to a side peak 1/f0 away).
% Change  : W (kernel length), B6 (bandwidth), rho (speckle correlation),
%           nTrial. A shorter kernel or narrower band moves the knee to
%           higher SNR; rho < 1 adds a high-SNR floor.
% Expect  : f0 = 5 MHz, B6 = 0.6, T = 2 us, rho = 1. From 15 to 30 dB both
%           estimators sit within about 1.3x of the Carter CRLB (0.41 ns
%           at 30 dB). At 10 dB a few % of trials hop to a side peak and the
%           std is about 10x the bound; at 5 dB about a third hop. Above
%           35 dB a floor near 0.3 ns appears (parabolic bias plus window-
%           edge decorrelation). For this Gaussian spectrum the Carter
%           integral and the WT closed form evaluated at B = B6 agree to
%           within 5 %.
rng(1);

%% Parameters
f0 = 5e6;  fs = 40e6;  B6 = 0.6;      % centre frequency, sampling rate (Hz), -6 dB fractional BW
N  = 512;  n0 = 200;   W = 80;        % record, kernel start (1-based), kernel length (samples)
K  = 24;                               % search range +-K lags
rho = 1;                               % speckle correlation between the records
snrdB  = -10:5:40;                     % SNR inside the -6 dB band (dB)
nTrial = 300;                          % trials per SNR
T = W / fs;                            % kernel duration (s)

%% Spectrum (one-sided analytic construction gives RF and IQ together)
f  = (0:N-1)' * fs / N;  pos = (f > 0) & (f < fs/2);   % positive-frequency bins
sf = B6 * f0 / 2 / sqrt(2*log(2));                       % spectral std of the amplitude response (Hz)
H  = exp(-(f - f0).^2 / (2*sf^2)) .* pos;                % amplitude response, positive bins only
Hn = H / sqrt(sum(H.^2) / N^2 * 2);                      % scale so the real part has unit power
nb = sum(pos & abs(f - f0) <= B6*f0/2);                  % bins inside the -6 dB band
fitPar = @(Rm, R0, Rp) 0.5 * (Rm - Rp) ./ (Rm - 2*R0 + Rp);
idxW = (n0:n0+W-1)';  lags = -K:K;
Yidx = bsxfun(@plus, idxW, lags);                        % W x (2K+1) window indices

sdPar = zeros(size(snrdB)); sdPha = sdPar; outFrac = sdPar; madPar = sdPar;
for is = 1:numel(snrdB)
    snr = 10^(snrdB(is)/10);
    sn  = sqrt(N^2 / (2 * nb * snr));                    % noise std per positive bin: in-band real-part noise power = 1/snr
    ePar = zeros(nTrial, 1); ePha = ePar;
    for t = 1:nTrial
        d  = rand;                                       % true shift, uniform in [0,1) samples
        S  = (randn(N,1) + 1j*randn(N,1)) / sqrt(2) .* Hn;
        U  = (randn(N,1) + 1j*randn(N,1)) / sqrt(2) .* Hn;
        Sy = (rho*S + sqrt(1 - rho^2)*U) .* exp(-1j*2*pi*f*d/fs);
        n1 = (randn(N,1) + 1j*randn(N,1)) / sqrt(2) * sn .* pos;
        n2 = (randn(N,1) + 1j*randn(N,1)) / sqrt(2) * sn .* pos;
        xa = ifft(2 * (S  + n1));  ya = ifft(2 * (Sy + n2));   % analytic signals; RF = real part
        x  = real(xa);  y = real(ya);
        xw = x(idxW);  Yw = y(Yidx);
        R  = (xw' * Yw) ./ sqrt(sum(xw.^2) * sum(Yw.^2, 1));   % NCC at all lags
        [~, im] = max(R(2:end-1));  im = im + 1;  k0 = lags(im);
        ePar(t) = k0 + fitPar(R(im-1), R(im), R(im+1)) - d;
        Cc = sum(ya(idxW + k0) .* conj(xa(idxW)));              % complex correlation at k0
        ePha(t) = k0 - angle(Cc) * fs / (2*pi*f0) - d;
    end
    ePar = ePar / fs;  ePha = ePha / fs;                  % seconds
    sdPar(is) = std(ePar);  sdPha(is) = std(ePha);
    madPar(is) = 1.4826 * median(abs(ePar - median(ePar)));
    outFrac(is) = mean(abs(ePar) > 0.5 / f0);
end

%% Bounds
% Carter CRLB: var >= [2T int (2 pi f)^2 g^2/(1-g^2) df]^-1, g = rho*S(f)/(S(f)+N0)
fg = linspace(1, fs/2, 4000);  dfg = fg(2) - fg(1);
Sg = exp(-(fg - f0).^2 / sf^2);  Sg = Sg / (sum(Sg) * dfg);     % unit-power signal PSD (one-sided)
crlb = zeros(size(snrdB)); wtb = crlb;
for is = 1:numel(snrdB)
    snr = 10^(snrdB(is)/10);  N0 = 1 / (snr * B6 * f0);         % white noise PSD: in-band noise power = 1/snr
    g = rho * Sg ./ (Sg + N0);
    crlb(is) = 1 / sqrt(2 * T * sum((2*pi*fg).^2 .* g.^2 ./ (1 - g.^2 + eps)) * dfg);
    wtb(is) = sqrt(3 / (2 * f0^3 * pi^2 * T * (B6^3 + 12*B6)) * ((1 + 1/snr)^2 / rho^2 - 1));
end

%% Print
fprintf('f0 = %.1f MHz, B6 = %.0f %%, T = %.2f us (%d smp), rho = %.3f, %d trials/SNR\n', f0/1e6, 100*B6, T*1e6, W, rho, nTrial);
fprintf(' SNR(dB)  sd_par(ns)  sd_phase(ns)  CRLB(ns)  WT(ns)  sd_par/CRLB  outliers\n');
for is = 1:numel(snrdB)
    fprintf('  %4d   %9.3f   %10.3f   %8.3f  %6.3f   %8.2f    %5.1f %%\n', snrdB(is), sdPar(is)*1e9, sdPha(is)*1e9, ...
        crlb(is)*1e9, wtb(is)*1e9, sdPar(is)/crlb(is), 100*outFrac(is));
end
fprintf('uniform over +-K lags would give sd = %.1f ns\n', 2*K/fs/sqrt(12)*1e9);

%% Figure
figure;
subplot(1, 2, 1);
semilogy(snrdB, sdPar*1e9, 'o-', snrdB, sdPha*1e9, 's--', snrdB, crlb*1e9, 'k-', snrdB, wtb*1e9, 'k:', snrdB, madPar*1e9, 'x');
xlabel('in-band SNR (dB)'); ylabel('std of delay error (ns)'); grid on;
legend('NCC+parab.', 'phase', 'CRLB', 'WT (B = B6)', '1.48 MAD', 'Location', 'northeast');
title('Jitter vs SNR');
subplot(1, 2, 2);
bar(snrdB, 100*outFrac); set(gca, 'XTick', -10:10:40); xlabel('SNR (dB)'); ylabel('outliers (%)'); grid on;
title('Outliers: |error| > 1/(2 f0)');
