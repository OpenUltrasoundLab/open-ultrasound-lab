%% m13_tikhonov_deconv_lcurve.m
% Tikhonov deconvolution of an RF line and four ways to choose lambda.
% Purpose : RF = pulse * reflectivity + noise (circulant model, so the SVD is
%           the DFT and sigma_i = |H(f_i)|). Tikhonov in Hansen's convention
%              x_lam = argmin ||A x - b||^2 + lam^2 ||x||^2,
%              filter factors f_i = sigma_i^2 / (sigma_i^2 + lam^2),
%           so lam is directly comparable to the singular values (|H| <= 1).
%           Compare the oracle lambda (needs the truth) with the discrepancy
%           principle, GCV and the L-curve corner, for white and for
%           correlated (AR(1), rho = 0.9) noise.
% Change  : snr_dB, the pulse bandwidth bw, the AR coefficient rho_ar.
% Expect  : white noise: GCV and L-curve land near the oracle, the
%           discrepancy principle (tau = 1) over-smooths a little.
%           AR(1) noise: GCV and the L-curve pick a tiny lambda and the error
%           explodes (they assume white noise); the discrepancy principle,
%           given the true noise norm, still works.
% Runtime : about 1 s.
clear; close all; rng(1);

%% Pulse and forward operator (circulant, applied by FFT)
N = 512; fs = 40e6; f0 = 5e6; bw = 0.6;        % samples, Hz, Hz, -6 dB fractional bandwidth
sf = bw*f0/2/sqrt(2*log(2));                   % spectral std of the Gaussian amplitude spectrum (Hz)
st = 1/(2*pi*sf);                              % envelope std (s)
t  = ((0:N-1)' - N/2)/fs;
h  = exp(-t.^2/(2*st^2)).*cos(2*pi*f0*t);
h  = circshift(h, -N/2);                       % pulse centred at sample 0 (zero-phase)
H  = fft(h); H = H/max(abs(H));                % normalise so sigma_max = 1
f  = (0:N-1)'*fs/N; f(f >= fs/2) = f(f >= fs/2) - fs;   % frequency axis (Hz)
Aop  = @(v) real(ifft(H.*fft(v)));

%% Reflectivity: a few interfaces plus weak diffuse scattering
x = 0.08*randn(N,1);
pos = [70 150 160 260 330 345 420]; amp = [1 -0.7 0.6 0.9 -0.5 0.5 0.8];
x(pos) = x(pos) + amp(:);
band = abs(H) >= 0.01;                         % band within 40 dB of the peak
X = fft(x); X(~band) = 0; xB = real(ifft(X));  % recoverable (band-limited) truth
fprintf('Recoverable band (|H| >= 0.01): %.2f to %.2f MHz, holds %.0f %% of the reflectivity energy\n', ...
    min(f(band & f > 0))/1e6, max(f(band & f > 0))/1e6, 100*norm(xB)^2/norm(x)^2);
b0 = Aop(x);

lams = logspace(-4, 0.5, 200);
snr_dB = 25; rho_ar = 0.9;
names = {'white', 'AR(1) rho=0.9'};
for nt = 1:2
    w = randn(N,1);
    if nt == 2, n = filter(1, [1 -rho_ar], w); else, n = w; end
    n = n*norm(b0)/norm(n)*10^(-snr_dB/20);    % set SNR
    sig = norm(n)/sqrt(N);                     % noise std (known for the discrepancy principle)
    b = b0 + n; B = fft(b);
    E = zeros(size(lams)); R = E; Xn = E; G = E;
    for k = 1:numel(lams)
        F  = abs(H).^2./(abs(H).^2 + lams(k)^2);           % filter factors
        xl = real(ifft(conj(H)./(abs(H).^2 + lams(k)^2).*B));
        r  = b - Aop(xl);
        E(k) = norm(xl - xB)/norm(xB);  R(k) = norm(r);  Xn(k) = norm(xl);
        G(k) = N*sum(r.^2)/(N - sum(F))^2;                 % GCV function
    end
    [~, io] = min(E); [~, ig] = min(G);
    [~, id] = min(abs(R - sqrt(N)*sig));                   % discrepancy, tau = 1
    lr = log(R); lx = log(Xn); tl = log(lams);             % L-curve curvature (log-log)
    d1r = gradient(lr, tl); d1x = gradient(lx, tl); d2r = gradient(d1r, tl); d2x = gradient(d1x, tl);
    kap = (d1r.*d2x - d2r.*d1x)./(d1r.^2 + d1x.^2).^1.5;
    [~, il] = max(kap);
    fprintf('\n%s noise, SNR %d dB (noise std %.4f)\n', names{nt}, snr_dB, sig);
    fprintf('  %-22s lambda = %.4f   rel. error = %.3f\n', 'oracle (needs truth)', lams(io), E(io));
    fprintf('  %-22s lambda = %.4f   rel. error = %.3f\n', 'discrepancy (tau=1)', lams(id), E(id));
    fprintf('  %-22s lambda = %.4f   rel. error = %.3f\n', 'GCV', lams(ig), E(ig));
    fprintf('  %-22s lambda = %.4f   rel. error = %.3f\n', 'L-curve corner', lams(il), E(il));
    if nt == 1
        Rw = R; Xw = Xn; Ew = E; Gw = G; iw = [io id ig il]; bw_ = b;
    end
end

%% Figures (white-noise case)
figure('Position', [50 50 1000 700]);
subplot(2,2,1);
loglog(Rw, Xw, 'k-'); hold on;
mk = {'o','s','d','^'}; cl = {'b','m','r',[0 0.6 0.4]};
for q = 1:4, loglog(Rw(iw(q)), Xw(iw(q)), mk{q}, 'Color', cl{q}, 'MarkerFaceColor', cl{q}); end
xlabel('residual norm ||b - A x_\lambda||'); ylabel('solution norm ||x_\lambda||');
title('L-curve (white noise)'); legend('L-curve','oracle','discrepancy','GCV','corner','Location','northeast');
subplot(2,2,2);
semilogx(lams, Ew, 'b-'); hold on;
semilogx(lams, Gw/max(Gw)*max(Ew), 'r--');
xlabel('\lambda (units of \sigma_{max})'); ylabel('relative error vs band-limited truth');
title('error (blue) and GCV (red, scaled)'); ylim([0 1.5]);
xl = real(ifft(conj(H)./(abs(H).^2 + lams(iw(1))^2).*fft(bw_)));
tt = (0:N-1)/fs*1e6;
subplot(2,2,3);
plot(tt, xB, 'Color', [0.6 0.6 0.6]); hold on; plot(tt, xl, 'b');
xlabel('time (\mus)'); ylabel('reflectivity'); title('band-limited truth (grey) and oracle Tikhonov (blue)');
subplot(2,2,4);
fp = f(1:N/2)/1e6; F = abs(H(1:N/2)).^2./(abs(H(1:N/2)).^2 + lams(iw(1))^2);
plot(fp, abs(H(1:N/2)), 'k--'); hold on; plot(fp, F, 'b');
xlabel('frequency (MHz)'); ylabel('gain'); title('|H(f)| (dashed) and filter factors at oracle \lambda');
