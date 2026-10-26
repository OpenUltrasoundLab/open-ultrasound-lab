%% m11_kernel_strain.m
% Purpose : Monte Carlo check of the kernel-length tradeoff under axial
%           strain (Demo 11.5). Point scatterers are imaged with a Gaussian
%           RF pulse, then compressed by a uniform strain eps about the
%           kernel centre and shifted by a random subsample delay d. NCC with
%           a cosine subsample fit estimates the delay at the kernel centre
%           for a range of kernel lengths T. The rms error is compared with
%           the Walker-Trahey bound evaluated with rho(eps, T), the
%           window-averaged RF correlation under strain.
% Change  : epsStrain (0 gives a monotone T^-1/2 curve), snrdB, f0, B6.
% Expect  : f0 = 5 MHz, B6 = 0.6, 20 dB in-band SNR, eps = 1 %: the bound
%           has its minimum at T = 1.5 us (2.09 ns = 1.6 um). The Monte Carlo
%           rms error has its minimum at T = 1-1.5 us (about 2.5 ns, 1.2x the
%           bound) and rises to 3.9 ns at 5 us as strain decorrelates the
%           kernel. Below 1 us side-peak hops (0.5-7 % of trials) make the
%           rms error 5-15x the bound while the robust spread (1.48 MAD)
%           still follows it: the threshold effect. With eps = 0 the error
%           falls as T^-1/2 (0.80 ns at 5 us, bound 0.81 ns).
rng(1);

%% Parameters
c  = 1540;  f0 = 5e6;  fs = 40e6;  B6 = 0.6;     % (m/s), (Hz), (Hz), -6 dB fractional BW
epsStrain = 0.01;                                 % axial strain (compression), 1 %
snrdB = 20;                                       % SNR inside the band (signal/noise power)
Tlist = [0.25 0.5 0.75 1 1.5 2 3 4 5] * 1e-6;     % kernel lengths (s)
nTrial = 200;                                     % speckle realizations
K = 12;                                           % search +-K lags
N = 640;  tc = 320;                               % record length and kernel centre (samples, 1-based)
dens = 40;                                        % scatterers per microsecond of record

%% Pulse and helpers
sigT = sqrt(2*log(2)) / (pi*B6*f0);               % Gaussian envelope std (s) for -6 dB BW = B6*f0
pulse = @(t) exp(-t.^2 / (2*sigT^2)) .* cos(2*pi*f0*t);
t = ((1:N)' - tc) / fs;                           % time axis, 0 at the kernel centre (s)
Ns = round(dens * N / fs * 1e6);                  % scatterers per record
fitCos = @(Rm, R0, Rp) atan((Rp - Rm) ./ (2*R0 .* sin(acos((Rm + Rp) ./ (2*R0))))) ./ acos((Rm + Rp) ./ (2*R0));
fitPar = @(Rm, R0, Rp) 0.5 * (Rm - Rp) ./ (Rm - 2*R0 + Rp);
lags = -K:K;
% band-limited noise: white noise filtered by the same pulse, so the coherence is flat across the band
hn = pulse((-60:60)' / fs);
err = zeros(nTrial, numel(Tlist));

for it = 1:nTrial
    tj = (rand(1, Ns) - 0.5) * N / fs;            % scatterer times (s), centred on the kernel
    aj = randn(1, Ns);                            % scattering amplitudes
    d  = rand / fs;                               % true delay at the kernel centre (s)
    x  = pulse(bsxfun(@minus, t, tj)) * aj';
    y  = pulse(bsxfun(@minus, t, tj*(1 - epsStrain) + d)) * aj';   % compressed about tc, then delayed
    sc = sqrt(mean(x.^2));
    n1 = conv(randn(N + 120, 1), hn, 'valid');  n1 = n1(1:N);
    n2 = conv(randn(N + 120, 1), hn, 'valid');  n2 = n2(1:N);
    gn = sc * 10^(-snrdB/20) / sqrt(mean([n1; n2].^2));
    x = x + gn*n1;  y = y + gn*n2;
    for iT = 1:numel(Tlist)
        W = round(Tlist(iT) * fs);
        idxW = (tc - floor(W/2) : tc - floor(W/2) + W - 1)';
        Yw = y(bsxfun(@plus, idxW, lags));  xw = x(idxW);
        R = (xw' * Yw) ./ sqrt(sum(xw.^2) * sum(Yw.^2, 1));
        [~, im] = max(R(2:end-1));  im = im + 1;
        dl = fitCos(R(im-1), R(im), R(im+1));
        if ~isfinite(dl) || abs(dl) > 1, dl = fitPar(R(im-1), R(im), R(im+1)); end
        err(it, iT) = (lags(im) + dl) / fs - d;   % kernel centre delay is d (strain is about tc)
    end
end

%% Bound: Walker-Trahey with rho(eps, T) = mean over the window of the RF ACF at lag eps*t
acf = @(tau) exp(-tau.^2 / (4*sigT^2)) .* cos(2*pi*f0*tau);
snr = 10^(snrdB/10);
Tf = linspace(0.1e-6, 6e-6, 300);  wt = zeros(size(Tf));
for i = 1:numel(Tf)
    tt = linspace(-Tf(i)/2, Tf(i)/2, 201);  rho = max(mean(acf(epsStrain * tt)), 1e-3);
    wt(i) = sqrt(3 / (2*f0^3*pi^2*Tf(i)*(B6^3 + 12*B6)) * ((1 + 1/snr)^2 / rho^2 - 1));
end
[wmin, iw] = min(wt);
rmsE = sqrt(mean(err.^2, 1));
hop  = mean(abs(err) > 0.5/f0, 1);
madE = 1.4826 * median(abs(err - median(err, 1)), 1);

%% Print
fprintf('f0 = %.1f MHz, B6 = %.0f %%, strain = %.1f %%, in-band SNR = %d dB, %d trials\n', f0/1e6, 100*B6, 100*epsStrain, snrdB, nTrial);
fprintf('bound minimum: T = %.2f us, sigma = %.2f ns (%.2f um)\n', Tf(iw)*1e6, wmin*1e9, wmin*c/2*1e6);
fprintf('  T(us)  periods  rms err(ns)  1.48 MAD(ns)  WT bound(ns)  hops\n');
for iT = 1:numel(Tlist)
    tt = linspace(-Tlist(iT)/2, Tlist(iT)/2, 201);  rho = mean(acf(epsStrain * tt));
    b = sqrt(3 / (2*f0^3*pi^2*Tlist(iT)*(B6^3 + 12*B6)) * ((1 + 1/snr)^2 / rho^2 - 1));
    fprintf('  %4.2f   %5.2f    %9.2f    %9.2f    %9.2f    %5.1f %%\n', Tlist(iT)*1e6, Tlist(iT)*f0, rmsE(iT)*1e9, madE(iT)*1e9, b*1e9, 100*hop(iT));
end

%% Figure
figure;
loglog(Tlist*1e6, rmsE*1e9, 'o-', Tlist*1e6, madE*1e9, 's--', Tf*1e6, wt*1e9, 'k-');
hold on; loglog([2 2]/f0*1e6, [0.1 1000], 'k:'); hold off;
xlabel('kernel length T (us)'); ylabel('delay error (ns)'); grid on; axis([0.1 6 0.3 300]);
legend('rms (Monte Carlo)', '1.48 MAD', 'WT bound with rho(eps,T)', 'T = 2/f0', 'Location', 'northeast');
title(sprintf('Kernel length vs strain: eps = %.1f %%, SNR = %d dB', 100*epsStrain, snrdB));
