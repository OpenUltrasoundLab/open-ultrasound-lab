%% m11_ncc_subsample_bias.m
% Purpose : bias of subsample delay estimators on band-limited RF speckle.
%           Two records share the same speckle; the second is delayed by a
%           fractional shift d (0..1 sample) with an exact FFT phase shift.
%           NCC is computed at integer lags, then refined with parabolic,
%           cosine, Gaussian and upsampled-sinc fits.
% Change  : f0, fs (fs/f0 sets the parabolic bias, roughly as (f0/fs)^2),
%           B6 (-6 dB fractional bandwidth), kernel length W.
% Expect  : noise-free, so the curves are bias only. At f0 = 5 MHz,
%           fs = 40 MHz, B6 = 0.6: max |bias| about 0.012 smp (parabolic),
%           0.025 (Gaussian on RF), 0.002 (cosine), < 0.001 (upsample x8).
%           Every curve is zero at d = 0, 0.5 and 1 (within the Monte Carlo
%           standard error) and odd about d = 0.5. The noise-free 'std over
%           speckle' (0.011 smp parabolic, 0.007 cosine, 0.0008 upsample)
%           shows that each kernel's interpolation error depends on its own
%           local spectrum, so the mean bias curve cannot simply be
%           subtracted as a correction.
rng(1);

%% Parameters
c   = 1540;          % speed of sound (m/s)
f0  = 5e6;           % centre frequency (Hz)
fs  = 40e6;          % sampling rate (Hz)
B6  = 0.6;           % -6 dB fractional bandwidth of the amplitude spectrum
N   = 512;           % record length (samples)
n0  = 200;           % kernel start (sample index, 1-based)
W   = 80;            % kernel length (samples) = 2 us
K   = 24;            % search range +-K lags
L   = 12;            % half-length of the Hann-tapered sinc interpolator (lags)
U   = 8;             % upsampling factor for the sinc interpolator
nReal = 400;         % speckle realizations averaged per shift
dGrid = linspace(0, 1, 21);   % true fractional shifts (samples)

%% Band-limited speckle generator (Gaussian amplitude spectrum)
f   = (0:N-1)' * fs / N;  f(f >= fs/2) = f(f >= fs/2) - fs;   % signed bin frequencies (Hz)
sf  = B6 * f0 / 2 / sqrt(2*log(2));                       % spectral std (Hz)
H   = exp(-(abs(f) - f0).^2 / (2*sf^2));                  % even amplitude response
H(N/2+1) = 0;                                             % no energy at Nyquist
snc = @(x) (sin(pi*x) + (x == 0)) ./ (pi*x + (x == 0));   % sinc without toolbox

%% Estimators applied to three NCC samples Rm, R0, Rp around the integer peak
fitPar = @(Rm, R0, Rp) 0.5 * (Rm - Rp) ./ (Rm - 2*R0 + Rp);
fitCos = @(Rm, R0, Rp) atan((Rp - Rm) ./ (2*R0 .* sin(acos((Rm + Rp) ./ (2*R0))))) ./ acos((Rm + Rp) ./ (2*R0));
fitGau = @(Rm, R0, Rp) 0.5 * (log(Rm) - log(Rp)) ./ (log(Rm) - 2*log(R0) + log(Rp));

lags = -K:K;
est  = zeros(numel(dGrid), nReal, 4);     % parabolic, cosine, Gaussian, upsampled sinc
idxW = (n0:n0+W-1)';  Yidx = bsxfun(@plus, idxW, lags);
for r = 1:nReal
    S  = fft(randn(N, 1)) .* H;           % speckle spectrum
    x  = real(ifft(S));
    for id = 1:numel(dGrid)
        d = dGrid(id);
        y = real(ifft(S .* exp(-1j*2*pi*f*d/fs)));   % exact fractional delay
        % NCC: fixed kernel of x, sliding window of y
        xw = x(idxW);  Yw = y(Yidx);                       % W x (2K+1) sliding windows
        R  = (xw' * Yw) ./ sqrt(sum(xw.^2) * sum(Yw.^2, 1));
        [~, im] = max(R(2:end-1));  im = im + 1;  k0 = lags(im);
        Rm = R(im-1); R0 = R(im); Rp = R(im+1);
        est(id, r, 1) = k0 + fitPar(Rm, R0, Rp);
        est(id, r, 2) = k0 + fitCos(Rm, R0, Rp);
        if min([Rm R0 Rp]) > 0, est(id, r, 3) = k0 + fitGau(Rm, R0, Rp); else, est(id, r, 3) = NaN; end
        % band-limited (sinc) interpolation of R on a 1/U grid, parabolic fit on the fine grid
        tf = k0 + (-U:U) / U;                              % fine lag grid (samples)
        u  = bsxfun(@minus, tf, lags');                    % lag offsets (numel(lags) x numel(tf))
        Rf = R * (snc(u) .* (0.5 + 0.5*cos(pi*u/(L+1))) .* (abs(u) <= L));
        [~, j] = max(Rf(2:end-1));  j = j + 1;
        est(id, r, 4) = tf(j) + fitPar(Rf(j-1), Rf(j), Rf(j+1)) / U;
    end
end

%% Bias = mean error over realizations; jitter from finite-window speckle
err  = bsxfun(@minus, est, dGrid');                       % (shift, realization, method)
bias = zeros(numel(dGrid), 4); jit = zeros(1, 4);   % jit: spread over speckle (finite kernel)
for m = 1:4
    e = err(:, :, m);
    for id = 1:numel(dGrid), v = e(id, :); bias(id, m) = mean(v(~isnan(v))); end
    v = e(:); v = v(~isnan(v)); jit(m) = std(v);
end
names = {'parabolic', 'cosine', 'Gaussian', 'upsample x8'};
fprintf('f0 = %.1f MHz, fs = %.1f MHz (fs/f0 = %.1f), B6 = %.0f %%, W = %d smp\n', f0/1e6, fs/1e6, fs/f0, 100*B6, W);
fprintf('1 sample = %.1f ns = %.2f um axial\n', 1e9/fs, c/(2*fs)*1e6);
for m = 1:4
    [mx, i] = max(abs(bias(:, m)));
    fprintf('%-12s max|bias| = %.4f smp (%.3f ns) at d = %.2f | std over speckle = %.4f smp\n', ...
        names{m}, mx, mx*1e9/fs, dGrid(i), jit(m));
end
fprintf('bias at d = 0 / 0.5 (parabolic): %.1e / %.1e smp (Monte Carlo standard error %.1e)\n', ...
    bias(1,1), bias(11,1), jit(1)/sqrt(nReal));

%% Figure
figure;
subplot(1, 2, 1);
plot(dGrid, bias(:,1), 'o-', dGrid, bias(:,2), 's-', dGrid, bias(:,3), '^-', dGrid, bias(:,4), 'd-');
hold on; plot([0 1], [0 0], 'k:'); hold off;
xlabel('true fractional shift d (samples)'); ylabel('bias (samples)');
title('Bias = mean estimate - true'); legend(names, 'Location', 'southoutside'); grid on;
subplot(1, 2, 2);
lagAx = (-24:24); tt = lagAx / fs; st = 1 / (2*pi*sf);
Rth = exp(-tt.^2 / (4*st^2)) .* cos(2*pi*f0*tt);          % expected RF autocorrelation
tfine = linspace(-24, 24, 2000) / fs;
plot(tfine*1e9, exp(-tfine.^2/(4*st^2)) .* cos(2*pi*f0*tfine), 'Color', [0.6 0.6 0.6]); hold on;
plot(tt*1e9, Rth, 'o'); hold off;
xlabel('lag (ns)'); ylabel('expected NCC'); title('Expected RF NCC'); xlim([-400 400]); grid on;
