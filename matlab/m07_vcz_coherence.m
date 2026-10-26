%% m07_vcz_coherence.m  --  van Cittert-Zernike: channel covariance and coherence vs lag
% Purpose : narrowband channel data from diffuse scatterers in the focal plane of a
%           64-element array (exact spherical geometry, focused transmit and receive).
%           Estimate the sample spatial covariance, show it as an image, and compare the
%           coherence-vs-lag curve with the VCZ triangle (autocorrelation of the transmit
%           aperture). Then add incoherent noise (lag >= 1 scaled by 1/(1+1/SNR)), a
%           near-field phase screen (two-way loss ~ exp(-D_phi(m))), compute SLSC and
%           the coherence factor, and check the Reed-Mallett-Brennan snapshot rule.
% Change  : Ntx (transmit aperture, elements), SNR_dB, sig_phi, ell_phi, M (snapshots).
% Expect  : clean C(1), C(16), C(32), C(48) = 0.99, 0.73, 0.48, 0.24 vs the paraxial
%           triangle 0.984, 0.750, 0.500, 0.250 (the +-6 mm scatterer window, a stand-in
%           for range gating, rounds the kink at lag 0; exact geometry at F/1.6 costs ~0.02);
%           0 dB noise halves lags >= 1; SLSC (lags 1-6) ~ 5.67 clean;
%           CF 0.649 (ratio of means; the paraxial model gives 2/3 + 1/(3N^2) = 0.667, the
%           exact geometry and scatterer window lower it) and 0.556 (mean of per-snapshot
%           ratios; 0.567 in the paraxial model): E[X/Y] differs from E[X]/E[Y];
%           Ntx = 32: C(16) = 0.51, C(32) = 0.01 (triangle 0.50, 0), the receive aperture
%           only samples the transmit autocorrelation;
%           RMB mean SINR loss ~ (K-N+2)/(K+1): 0.50 at K = 2N-3.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);

%% 1. Geometry
c0 = 1540;  f0 = 5e6;  lam = c0/f0;  k = 2*pi/lam;   % m/s, Hz, m, rad/m
N  = 64;  pitch = 0.3e-3;  z = 30e-3;                 % elements, m, focal depth (m)
xe = ((0:N-1) - (N-1)/2) * pitch;                     % element x positions (m), 1 x N
Ntx = 64;                                             % active transmit elements (centred); try 32
tx  = abs(((0:N-1) - (N-1)/2)) < Ntx/2;               % transmit apodization (uniform), 1 x N
xi = (-6e-3 : 20e-6 : 6e-3);                          % scatterer grid across the focal plane (m)
Nxi = numel(xi);
M  = 2000;                                            % independent speckle realizations

% one-way propagation phasors element -> scatterer, minus the focusing delay to (0,z)
r  = sqrt((xe.' - xi).^2 + z^2);                      % N x Nxi distances (m)
r0 = sqrt(xe.'.^2 + z^2);                             % N x 1 distances to the focus
G  = exp(-1i*k*(r - r0));                             % focused one-way kernel, N x Nxi

Pt = (tx * G);                                        % transmit field at each scatterer, 1 x Nxi
B  = G .* Pt;                                         % receive kernel weighted by transmit field
in6 = xi(abs(Pt) >= 0.5*max(abs(Pt)));               % -6 dB (amplitude) extent of the transmit beam
fprintf('transmit -6 dB beam width %.2f mm (uniform aperture theory 1.21 lambda F# = %.2f mm)\n', ...
    1e3*(max(in6) - min(in6)), 1e3*1.207*lam*z/(Ntx*pitch));

%% 2. Channel data, sample covariance, coherence vs lag
A  = (randn(Nxi, M) + 1i*randn(Nxi, M)) / sqrt(2);   % delta-correlated complex reflectivity
S0 = B * A;                                           % N x M focused channel samples
Ps = mean(abs(S0(:)).^2);                             % signal power per channel
cohfun = @(S) (S * S') / size(S, 2);                  % sample covariance R = (1/M) S S^H
nrm    = @(R) R ./ sqrt(real(diag(R)) * real(diag(R)).');
lagavg = @(C) arrayfun(@(m) mean(real(diag(C, m))), 0:size(C,1)-1);
tri = max(0, 1 - (0:N-1)/Ntx);                        % VCZ prediction: transmit autocorrelation

R0 = cohfun(S0);  C0 = nrm(R0);  coh0 = lagavg(C0);
fprintf('clean: C(1) %.3f  C(16) %.3f  C(32) %.3f  C(48) %.3f   (triangle %.3f %.3f %.3f %.3f)\n', ...
    coh0([2 17 33 49]), tri([2 17 33 49]));
disp('top-left 4x4 corner of the normalized sample covariance:'); disp(round(100*C0(1:4,1:4))/100);

%% 3. Incoherent noise: lag 0 stays 1, lags >= 1 scale by 1/(1+1/SNR)
SNR_dB = 0;  snr = 10^(SNR_dB/10);
Wn = sqrt(Ps/snr) * (randn(N, M) + 1i*randn(N, M)) / sqrt(2);
cohN = lagavg(nrm(cohfun(S0 + Wn)));
fprintf('SNR %g dB: C(1) %.3f, ratio to clean %.3f (theory 1/(1+1/SNR) = %.3f)\n', ...
    SNR_dB, cohN(2), cohN(2)/coh0(2), 1/(1 + 1/snr));

%% 4. Near-field phase screen on transmit and receive
sig_phi = 0.7;  ell_phi = 4;                          % rms phase (rad), correlation length (elements)
kk = -15:15;  g = exp(-kk.^2 / ell_phi^2);            % filter -> corr(m) = exp(-m^2/(2 ell^2))
w  = conv(randn(1, N + 30), g, 'same');  phi = w(16:15+N);
phi = phi - mean(phi);  phi = sig_phi * phi / sqrt(mean(phi.^2));
Gab = G .* exp(1i*phi.');                             % screen at the aperture, both ways
Bab = Gab .* (tx * Gab);
cohA = lagavg(nrm(cohfun(Bab * A)));
avgph = arrayfun(@(m) abs(mean(exp(1i*(phi(1:N-m) - phi(1+m:N)))))^2, 0:N-1);
Dphi  = arrayfun(@(m) mean((phi(1:N-m) - phi(1+m:N)).^2), 0:N-1);
fprintf('screen %.1f rad, ell %d el.: lag  measured  tri*|<e^jdphi>|^2  tri*exp(-D)\n', sig_phi, ell_phi);
for m = [1 2 4 8 16]
    fprintf('                          %2d    %.3f     %.3f              %.3f\n', ...
        m, cohA(m+1), tri(m+1)*avgph(m+1), tri(m+1)*exp(-Dphi(m+1)));
end

%% 5. SLSC (lags 1..6) and coherence factor
Mlag = 6;
slsc = @(coh) sum(coh(2:Mlag+1));
cf_ratio = @(S) mean(abs(sum(S, 1)).^2) / mean(N * sum(abs(S).^2, 1));   % ratio of means
cf_mean  = @(S) mean(abs(sum(S, 1)).^2 ./ (N * sum(abs(S).^2, 1)));      % mean of ratios
fprintf('SLSC: clean %.2f (ideal %.2f), 0 dB noise %.2f, screen %.2f\n', ...
    slsc(coh0), sum(tri(2:Mlag+1)), slsc(cohN), slsc(cohA));
fprintf('CF clean: ratio of means %.3f (paraxial theory %.3f), mean of ratios %.3f\n', ...
    cf_ratio(S0), sum((N - abs(-(N-1):(N-1))) .* tri(abs(-(N-1):(N-1)) + 1)) / N^2, cf_mean(S0));
fprintf('CF 0 dB noise %.3f, CF screen %.3f\n', cf_ratio(S0 + Wn), cf_ratio(Bab * A));

%% 6. Reed-Mallett-Brennan: sample-matrix MVDR needs K ~ 2N snapshots for 3 dB
Nr = 32;  Rt = toeplitz(1 - (0:Nr-1)/Nr) + 0.01*eye(Nr);   % speckle (VCZ) + 20 dB noise
Lr = chol(Rt, 'lower');  s = ones(Nr, 1);  opt = real(s' * (Rt \ s));
for K = [Nr, 2*Nr-3, 4*Nr]
    loss = zeros(1, 500);
    for t = 1:500
        X = Lr * (randn(Nr, K) + 1i*randn(Nr, K)) / sqrt(2);
        wv = (X*X'/K) \ s;                                    % sample-matrix-inverse weights
        loss(t) = abs(wv'*s)^2 / real(wv'*Rt*wv) / opt;       % SINR relative to optimum
    end
    fprintf('RMB N = %d, K = %3d: mean SINR loss %.3f (theory (K-N+2)/(K+1) = %.3f, %.2f dB)\n', ...
        Nr, K, mean(loss), (K-Nr+2)/(K+1), 10*log10((K-Nr+2)/(K+1)));
end

%% 7. Figures
figure('Position', [100 100 1050 420]);
subplot(1,2,1);
imagesc(1:N, 1:N, abs(C0)); axis image; colormap(gray); colorbar;
xlabel('element j'); ylabel('element i'); title('|R_{ij}| / sqrt(R_{ii} R_{jj}), clean speckle');
subplot(1,2,2);
plot(0:N-1, tri, 'k--', 'LineWidth', 1.5); hold on;
plot(0:N-1, coh0, 'b-', 'LineWidth', 1.5);
plot(0:N-1, cohN, 'r-', 'LineWidth', 1.2);
plot(0:N-1, cohA, 'Color', [0 0.6 0.4], 'LineWidth', 1.2);
plot(0:N-1, tri .* avgph, ':', 'Color', [0 0.6 0.4], 'LineWidth', 1.5);
grid on; xlim([0 N]); ylim([-0.05 1.05]);
xlabel('element lag m'); ylabel('coherence C(m)');
legend('VCZ triangle', 'clean', sprintf('SNR %g dB', SNR_dB), 'phase screen', 'screen theory', 'Location', 'northeast');
title('spatial coherence vs lag');
