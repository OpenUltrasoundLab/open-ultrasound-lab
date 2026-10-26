%% m27_multitaper_dpss.m  -  Slepian (DPSS) tapers from eig, and Thomson's multitaper spectrum
% Purpose : (1) build the K = 2NW - 1 discrete prolate spheroidal sequences of
%           length N without a toolbox, as eigenvectors of the symmetric
%           tridiagonal matrix that commutes with the sinc (prolate) kernel;
%           (2) check them against the sinc-kernel eigenproblem: the
%           concentration lambda_k (fraction of taper energy inside |f| <= W)
%           equals the sinc-kernel eigenvalue;
%           (3) estimate the spectrum of a short speckle gate with K tapers and
%           verify that the variance falls as 1/K (nu = 2K dof);
%           (4) compare variance x bandwidth with Hann-Welch at 50 % overlap.
% Change  : N, NW (half bandwidth in bins), K.
% Expect  : NW = 4, N = 128: lambda = 1.0000 ... 0.9994 0.9926 0.9369, then
%           lambda_7 = 0.70 (leaky); orthonormality error ~1e-15;
%           var(S_mt/S) = 0.17, a little above 1/K = 0.143 because the speckle
%           spectrum is not flat across 2W = 2NW*fs/N = 2.5 MHz; the same
%           smoothing biases the mean low (0.91) on this short gate. Hann-Welch
%           with the same spectral-window width gives about the same variance.
clear; close all; rng(1);

%% Parameters
N  = 128;            % gate length (samples); 3.2 us at 40 MHz = 2.46 mm of depth
NW = 4;              % time-half-bandwidth product
W  = NW/N;           % half bandwidth (cycles/sample)
K  = 2*NW - 1;       % number of well-concentrated tapers
fs = 40;             % MHz

%% 1. DPSS from the tridiagonal matrix (Slepian 1978; Percival & Walden 1993, ch. 8)
n = (0:N-1)';
d = ((N - 1)/2 - n).^2*cos(2*pi*W);          % diagonal
e = n(2:end).*(N - n(2:end))/2;              % off-diagonal
T = diag(d) + diag(e, 1) + diag(e, -1);
[V, D] = eig(T);  [~, idx] = sort(diag(D), 'descend');
V = V(:, idx(1:K+1));                        % keep one extra taper to show the leakage cliff
for k = 1:K+1                                % sign convention: even tapers sum > 0, odd tapers start > 0 slope
  if mod(k-1, 2) == 0, s = sign(sum(V(:, k))); else, s = sign(sum((n - (N-1)/2).*V(:, k))); end
  V(:, k) = s*V(:, k);
end
fprintf('Orthonormality: max |V''V - I| = %.1e\n', max(max(abs(V'*V - eye(K+1)))));

%% 2. Concentration from the sinc kernel  A(m,n) = sin(2 pi W (m-n)) / (pi (m-n))
[mm, nn] = meshgrid(n, n);  dmn = mm - nn;
A = sin(2*pi*W*dmn)./(pi*dmn + (dmn == 0));  A(dmn == 0) = 2*W;
lam = diag(V'*A*V);                          % Rayleigh quotient = energy fraction in |f| <= W
lamA = sort(eig((A + A')/2), 'descend');
fprintf('\n k   lambda_k (from tapers)   eig of sinc kernel   out-of-band energy (dB)\n');
for k = 1:K+1
  fprintf('%2d     %.6f             %.6f            %7.1f\n', k-1, lam(k), lamA(k), 10*log10(max(1 - lam(k), 1e-300)));
end
fprintf('Number of eigenvalues > 0.5: %d (Shannon number 2NW = %d)\n', sum(lamA > 0.5), 2*NW);

%% 3. Multitaper spectrum of a speckle gate, Monte Carlo variance
f0 = 5;  bw = 0.6;  sigf = bw*f0/(2*sqrt(2*log(2)));
Nf = 1024;  ff = (0:Nf-1)'/Nf*fs;  ff(ff >= fs/2) = ff(ff >= fs/2) - fs;
H = exp(-(abs(ff) - f0).^2/(2*sigf^2));
Strue = @(f) exp(-(abs(f) - f0).^2/sigf^2);
ntr = 400;  Np = 512;  fP = (0:Np-1)'/Np*fs;  kb = abs(fP - f0) <= 1.2;
Rmt = zeros(sum(kb), ntr);  Rhw = Rmt;  Rp = Rmt;
w = 0.5 - 0.5*cos(2*pi*(n + 0.5)/N);  w = w/sqrt(sum(w.^2));
M = 32;  S = 16;  st = 1:S:(N - M + 1);  wh = 0.5 - 0.5*cos(2*pi*((0:M-1)' + 0.5)/M);  wh = wh/sqrt(sum(wh.^2));
for t = 1:ntr
  x = real(ifft(fft(randn(Nf, 1)).*H));  x = x(1:N);
  Smt = mean(abs(fft(V(:, 1:K).*repmat(x, 1, K), Np)).^2, 2);
  Sp  = abs(fft(x.*w, Np)).^2;
  Sw  = zeros(Np, 1);
  for b = 1:numel(st), Sw = Sw + abs(fft(x(st(b):st(b)+M-1).*wh, Np)).^2; end
  Sw = Sw/numel(st);
  Rmt(:, t) = Smt(kb)./Strue(fP(kb));  Rp(:, t) = Sp(kb)./Strue(fP(kb));  Rhw(:, t) = Sw(kb)./Strue(fP(kb));
end
fprintf('\nSpeckle gate, N = %d (%.2f mm), bins within f0 +/- 1.2 MHz, %d trials:\n', N, N/fs*1.54/2, ntr);
fprintf('  Hann periodogram : var(S^/S) = %.3f, mean = %.3f\n', var(Rp(:)), mean(Rp(:)));
fprintf('  Multitaper K = %d: var(S^/S) = %.3f (1/K = %.3f), mean = %.3f\n', K, var(Rmt(:)), 1/K, mean(Rmt(:)));
fprintf('  Hann-Welch M = %d, 50%% overlap (%d segments): var = %.3f, mean = %.3f\n', M, numel(st), var(Rhw(:)), mean(Rhw(:)));
fprintf('  Resolution: multitaper 2W = %.3f MHz; Welch bin fs/M = %.3f MHz\n', 2*W*fs, fs/M);

%% 4. Bandwidth (int Wk)^2/int Wk^2 of the averaged spectral windows (in MHz)
L = 32*N;  fL = (0:L-1)'/L*fs;
Wmt = mean(abs(fft(V(:, 1:K), L)).^2, 2);  Wmt = Wmt/sum(Wmt);
Whw = abs(fft(wh, L)).^2;  Whw = Whw/sum(Whw);
Bmt = 1/sum(Wmt.^2)*fs/L;  Bhw = 1/sum(Whw.^2)*fs/L;
fprintf('  Spectral-window bandwidth: multitaper %.2f MHz, Welch segment %.2f MHz\n', Bmt, Bhw);
fprintf('  variance x bandwidth x T: multitaper %.2f, Welch %.2f (T = %.2f us)\n', ...
        var(Rmt(:))*Bmt*N/fs, var(Rhw(:))*Bhw*N/fs, N/fs);

%% Figures
figure('Position', [100 100 950 380]);
subplot(1, 2, 1);
plot(n, V(:, 1:4), 'LineWidth', 1.3); grid on; xlim([0 N-1]);
xlabel('sample n'); ylabel('taper value');
legend('v_0', 'v_1', 'v_2', 'v_3', 'Location', 'south');
title(sprintf('DPSS tapers, N = %d, NW = %g', N, NW));
subplot(1, 2, 2);
fW = fL;  fW(fW >= fs/2) = fW(fW >= fs/2) - fs;  [fW, o] = sort(fW);
Wk = abs(fft(V, L)).^2;  Wk = Wk(o, :);
plot(fW, 10*log10(Wk(:, 1)/max(Wk(:, 1))), 'b', 'LineWidth', 1.3); hold on;
plot(fW, 10*log10(Wk(:, K)/max(Wk(:, K))), 'm', 'LineWidth', 1.1);
plot(fW, 10*log10(Wk(:, K+1)/max(Wk(:, K+1))), 'Color', [0.6 0.6 0.6]);
plot([W W]*fs, [-120 5], 'k--', [-W -W]*fs, [-120 5], 'k--');
xlim([-6 6]); ylim([-120 5]); grid on;
xlabel('frequency offset (MHz)'); ylabel('|V_k(f)|^2 (dB)');
legend('k = 0', sprintf('k = %d', K-1), sprintf('k = %d (leaky)', K), '\pm W', 'Location', 'south');
title('Spectral windows: concentration in |f| \leq W');
