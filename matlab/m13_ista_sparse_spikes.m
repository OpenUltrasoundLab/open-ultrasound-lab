%% m13_ista_sparse_spikes.m
% Sparse-spike deconvolution: L2 (Tikhonov) vs L1 (ISTA and FISTA).
% Purpose : Recover point reflectors from a band-limited RF line by
%              min_x 0.5 ||A x - b||^2 + lam ||x||_1
%           with ISTA  x <- soft(x - A'(A x - b)/Lip, lam/Lip)
%           and FISTA (Beck & Teboulle 2009, Nesterov momentum), Lip = sigma_max^2 = 1.
%           Check the worst-case bounds  F(x_k) - F* <= Lip ||x0 - x*||^2 / (2k)    (ISTA)
%                                         F(x_k) - F* <= 2 Lip ||x0 - x*||^2 / (k+1)^2 (FISTA).
%           lam starts from the universal threshold sigma_n ||h|| sqrt(2 ln N).
%           Debias: least squares on the recovered support removes the L1 shrinkage.
% Change  : spike list, snr_dB, lam_factor, number of iterations.
% Expect  : FISTA reaches a 1e-6 relative gap in a few hundred iterations,
%           ISTA does not within 1000; L1 recovers the isolated spikes and the
%           opposite-sign pair 6 samples apart; the same-sign pair 4 samples apart
%           (0.1 us = half a period at 5 MHz, lambda/4 in depth) is lost: its two
%           echoes nearly cancel in band (1 - 0.9 = 0.1 at f0), so no prior
%           can recover it and L1 returns wrong-sign spikes at the wrong positions.
% Runtime : about 2 s.
clear; close all; rng(2);

N = 256; fs = 40e6; f0 = 5e6; bw = 0.6;
sf = bw*f0/2/sqrt(2*log(2)); st = 1/(2*pi*sf);
t = ((0:N-1)' - N/2)/fs;
h = circshift(exp(-t.^2/(2*st^2)).*cos(2*pi*f0*t), -N/2);
H = fft(h); H = H/max(abs(H)); h = real(ifft(H));
Aop  = @(v) real(ifft(H.*fft(v)));
Atop = @(v) real(ifft(conj(H).*fft(v)));
soft = @(v, th) sign(v).*max(abs(v) - th, 0);
fprintf('Pulse: -6 dB envelope duration %.0f ns = %.1f samples; ||h|| = %.3f\n', ...
    2*sqrt(2*log(2))*st*1e9, 2*sqrt(2*log(2))*st*fs, norm(h));

x = zeros(N,1);
pos = [30 70 76 120 170 174 215]; amp = [1 -0.8 0.7 0.6 -1 -0.9 0.5];
x(pos + 1) = amp(:);                             % 0-based sample positions
b0 = Aop(x); n = randn(N,1); snr_dB = 25;
n = n*norm(b0)/norm(n)*10^(-snr_dB/20); sig = norm(n)/sqrt(N); b = b0 + n;
lamU = sig*norm(h)*sqrt(2*log(N));
lam = 1.0*lamU;
Fobj = @(v) 0.5*sum((Aop(v) - b).^2) + lam*sum(abs(v));
fprintf('noise std %.4f, universal lambda %.4f, using lambda = %.4f\n', sig, lamU, lam);

%% Reference minimizer (long FISTA run)
xs = zeros(N,1); y = xs; tk = 1;
for k = 1:20000
    xn = soft(y - Atop(Aop(y) - b), lam); tn = (1 + sqrt(1 + 4*tk^2))/2;
    y = xn + ((tk - 1)/tn)*(xn - xs); xs = xn; tk = tn;
end
Fs = Fobj(xs); R2 = sum(xs.^2);

%% ISTA vs FISTA from x0 = 0
K = 1000; gi = zeros(K,1); gf = gi;
xi = zeros(N,1); xf = zeros(N,1); y = xf; tk = 1;
for k = 1:K
    xi = soft(xi - Atop(Aop(xi) - b), lam); gi(k) = Fobj(xi) - Fs;
    xn = soft(y - Atop(Aop(y) - b), lam); tn = (1 + sqrt(1 + 4*tk^2))/2;
    y = xn + ((tk - 1)/tn)*(xn - xf); xf = xn; tk = tn; gf(k) = Fobj(xf) - Fs;
end
kk = (1:K)'; bI = R2./(2*kk); bF = 2*R2./(kk + 1).^2;
fprintf('relative gap (F-F*)/F* at k = 10, 100, 1000:  ISTA %.1e %.1e %.1e   FISTA %.1e %.1e %.1e\n', ...
    gi(10)/Fs, gi(100)/Fs, gi(1000)/Fs, gf(10)/Fs, gf(100)/Fs, gf(1000)/Fs);
fprintf('worst-case bounds hold: ISTA %d, FISTA %d\n', all(gi <= bI + 1e-12), all(gf <= bF + 1e-12));
k6i = find(gi/Fs < 1e-6, 1); k6f = find(gf/Fs < 1e-6, 1);
if isempty(k6i), k6i = NaN; end
fprintf('iterations to a 1e-6 relative gap: ISTA %g, FISTA %d\n', k6i, k6f);

%% Support, debiasing, and the L2 comparison
supp = abs(xs) > 0.05*max(abs(xs));
xd = zeros(N,1); idx = find(supp);
Ad = zeros(N, numel(idx));
for q = 1:numel(idx), e = zeros(N,1); e(idx(q)) = 1; Ad(:,q) = Aop(e); end
xd(idx) = Ad\b;                                  % least squares on the support
fprintf('L1 support (|x| > 5%% of max): %d spikes at samples %s\n', numel(idx), mat2str(idx' - 1));
fprintf('L1 amplitudes on that support: %s\n', mat2str(round(xs(idx)'*100)/100));
fprintf('true spikes at samples %s\n', mat2str(pos));
fprintf('L1 amplitudes at true positions:       %s\n', mat2str(round(xs(pos + 1)'*100)/100));
fprintf('debiased amplitudes at true positions: %s\n', mat2str(round(xd(pos + 1)'*100)/100));
lams = logspace(-4, 0.5, 200); Et = zeros(size(lams)); B = fft(b);
for q = 1:numel(lams)
    Et(q) = norm(real(ifft(conj(H)./(abs(H).^2 + lams(q)^2).*B)) - x)/norm(x);
end
[et, qt] = min(Et); xt = real(ifft(conj(H)./(abs(H).^2 + lams(qt)^2).*B));
fprintf('relative error: Tikhonov (oracle lambda %.4f) %.3f | L1 %.3f | L1 debiased %.3f\n', ...
    lams(qt), et, norm(xs - x)/norm(x), norm(xd - x)/norm(x));
% (Indices printed 0-based to match the browser demo; MATLAB arrays are 1-based,
%  so x(pos+1) is the spike listed at sample pos.)

%% Figures
tt = (0:N-1)/fs*1e6;
figure('Position', [50 50 1000 650]);
subplot(2,2,1); plot(tt, b, 'Color', [0.8 0.5 0]); xlabel('time (\mus)'); title('RF data b');
subplot(2,2,2); stem(tt(pos + 1), amp, 'k', 'filled'); hold on; plot(tt, xt, 'm');
xlabel('time (\mus)'); title('truth (stems) and Tikhonov (oracle \lambda)');
subplot(2,2,3); stem(tt(pos + 1), amp, 'k', 'filled'); hold on; stem(tt(idx), xs(idx), 'b');
xlabel('time (\mus)'); title('truth (black) and L1 / FISTA (blue)');
subplot(2,2,4); loglog(kk, gi, 'r', kk, gf, 'b', kk, bI, 'r:', kk, bF, 'b:');
xlabel('iteration k'); ylabel('F(x_k) - F^*'); title('ISTA (red) vs FISTA (blue), dotted: bounds');
