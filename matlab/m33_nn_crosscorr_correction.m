%% m33_nn_crosscorr_correction.m
% Aberration law from speckle: nearest-neighbour cross-correlation (Flax and
% O'Donnell), its random-walk drift, multi-lag least squares, and one
% transmit iteration.
% Purpose : 128 elements, pitch 0.3 mm, 5 MHz, 60 % bandwidth, focus 40 mm.
%           A Gaussian delay screen (rms sig_tau, correlation length l) sits on
%           the array, on transmit and on receive. Speckle IQ data are
%           synthesized per element after the geometric focusing delays:
%             s_i(t) = sum_s b_s a_s env(t - t_is) exp(-j 2 pi f0 t_is),
%           t_is = residual geometric delay of scatterer s after dynamic receive
%           focusing (lateral offset only) + its two-way time + tau_i, and b_s the
%           (aberrated or corrected) transmit beam at the scatterer.
%           NN estimate : d_i = angle(sum_t s_i conj(s_i+1)) / (2 pi f0),
%                         tau_hat = cumulative sum of d_i (a random walk).
%           LS estimate : all pairs up to lag L, weighted by |rho_ij|^2/(1-|rho_ij|^2),
%                         solved for tau with sum(tau) = 0.
%           Mean and linear ramp are removed from truth and estimates before
%           comparison (a ramp only steers the beam).
%           Prediction for the detrended NN error: sigma_eps sqrt(N/15)
%           (Brownian motion minus its least-squares line).
% Change  : sig_tau, ell, Lmax, R, Ns.
% Expect  : noiseless speckle: sub-ns pair jitter, but a smooth bias of a few
%           ns from the blurred (aberrated) transmit focus; LS cannot remove a
%           bias, a second pass with the corrected transmit halves it.
%           0 dB channel SNR: jitter dominates, the NN error drifts as a random
%           walk (sigma_eps sqrt(N/15) after detrending is the bound for
%           independent pair errors; element noise enters adjacent pairs with
%           opposite signs and partly cancels) and LS over 8 lags cuts it.
% Runtime : about 12 s in Octave.
clear; close all; rng(2);

c = 1540; f0 = 5e6; lam = c/f0; k = 2*pi/lam;
N = 128; p = 0.3e-3; xe = ((0:N-1) - (N-1)/2)*p; D = N*p; zf = 40e-3;
sig_tau = 40e-9; ell = 4e-3; Lmax = 8; R = 16; Ns = 400;
snr_list = [Inf 0];                                          % channel SNR (dB): noiseless, and 0 dB per element
st = 1/(2*pi*0.6*f0/2.355);                                 % envelope std for 60 % (-6 dB) bandwidth
t = (-2e-6:0.05e-6:2e-6)';                                  % fast time around the focus (s), IQ at 20 MHz

dxs = p/4; xs = -60e-3:dxs:60e-3; g = exp(-2*(((-3*ell):dxs:(3*ell))/ell).^2); g = g/sqrt(sum(g.^2));
idx = round((xe - xs(1))/dxs) + 1;
detr = @(v) v - [ones(N, 1) xe(:)]*([ones(N, 1) xe(:)]\v(:));   % remove mean and linear ramp

% LS design matrix for all pairs with lag 1..Lmax (+ one row for sum(tau) = 0)
pairs = []; for L = 1:Lmax, pairs = [pairs; (1:N-L)' ((1:N-L) + L)']; end
A = zeros(size(pairs, 1), N); for q = 1:size(pairs, 1), A(q, pairs(q, 1)) = -1; A(q, pairs(q, 2)) = 1; end

TAU = zeros(N, R); XS = zeros(Ns, R); ZS = XS; AS = XS;      % the same screens and scatterers for every SNR
for r = 1:R
    tmp = sig_tau*conv(randn(1, numel(xs)), g, 'same'); TAU(:, r) = tmp(idx)';
    XS(:, r) = (rand(Ns, 1) - 0.5)*6e-3; ZS(:, r) = zf + (rand(Ns, 1) - 0.5)*4e-3;
    AS(:, r) = (randn(Ns, 1) + 1j*randn(Ns, 1))/sqrt(2);
end
for isnr = 1:numel(snr_list)
snr = snr_list(isnr);
drift = zeros(N, 1); res = zeros(R, 2, 3); eps_nn = zeros(R, 2); S_before = zeros(R, 1); S_after = zeros(R, 2, 2);
drift = zeros(N, 1);
for r = 1:R
    tau = TAU(:, r); phi = 2*pi*f0*tau; xsc = XS(:, r); zsc = ZS(:, r); a_s = AS(:, r);
    tx_law = zeros(N, 1);                                    % pass 1: transmit with the geometric law only
    S_before(r) = abs(mean(exp(1j*phi)))^2;
    for pass = 1:2
        b = exp(-1j*k*xsc*xe/zf)*exp(1j*(phi - tx_law))/N;   % transmit beam at each scatterer (Fraunhofer)
        b = b.*exp(-((zsc - zf)/2e-3).^2);                  % soft depth weighting of the kernel
        S = zeros(numel(t), N);
        for i = 1:N
            tis = (sqrt((xsc - xe(i)).^2 + zsc.^2) - sqrt(xe(i)^2 + zsc.^2) + 2*(zsc - zf))/c + tau(i);   % dynamic receive focus
            S(:, i) = (exp(-(t - tis').^2/(2*st^2)).*exp(-1j*2*pi*f0*tis'))*(b.*a_s);
        end
        if isfinite(snr), S = S + sqrt(mean(abs(S(:)).^2)/10^(snr/10))*(randn(size(S)) + 1j*randn(size(S)))/sqrt(2); end
        % NN estimate
        d = zeros(N-1, 1); rho1 = d;
        for i = 1:N-1
            q = sum(S(:, i).*conj(S(:, i+1))); d(i) = angle(q)/(2*pi*f0);    % estimate of tau_{i+1} - tau_i
            rho1(i) = abs(q)/sqrt(sum(abs(S(:, i)).^2)*sum(abs(S(:, i+1)).^2));
        end
        eps_nn(r, pass) = std(d - diff(tau));
        tnn = [0; cumsum(d)];
        % LS multi-lag estimate
        dd = zeros(size(pairs, 1), 1); w = dd;
        for q = 1:size(pairs, 1)
            i = pairs(q, 1); j = pairs(q, 2); cc = sum(S(:, i).*conj(S(:, j)));
            dd(q) = angle(cc)/(2*pi*f0); rh = min(abs(cc)/sqrt(sum(abs(S(:, i)).^2)*sum(abs(S(:, j)).^2)), 0.999);
            w(q) = rh^2/(1 - rh^2);
        end
        W = sqrt(w); tls = [A.*W; ones(1, N)]\[dd.*W; 0];
        et = detr(tau); e1 = detr(tnn) - et; e2 = detr(tls) - et;
        res(r, pass, :) = [std(et) std(e1) std(e2)];
        S_after(r, pass, :) = [abs(mean(exp(1j*2*pi*f0*e1)))^2 abs(mean(exp(1j*2*pi*f0*e2)))^2];
        if pass == 1, drift = drift + e1.^2; mr1(r) = mean(rho1); else, mr2(r) = mean(rho1); end
        tx_law = 2*pi*f0*tls;                                % pass 2: transmit with the LS estimate
    end
end
fprintf('\n===== channel SNR = %g dB =====\n', snr);
fprintf('Screen: rms %.0f ns, l = %.0f mm; %d elements; %d realizations\n', sig_tau*1e9, ell*1e3, N, R);
fprintf('rms of the detrended true screen: %.1f ns; Strehl of the raw screen %.3f\n', mean(res(:, 1, 1))*1e9, mean(S_before));
for pass = 1:2
    if pass == 1, lab = 'pass 1 (aberrated transmit)'; mr = mean(mr1); else, lab = 'pass 2 (transmit corrected by pass-1 LS)'; mr = mean(mr2); end
    fprintf('%s: mean adjacent-element |rho| = %.3f\n', lab, mr);
    fprintf('  NN pair jitter sigma_eps = %.2f ns\n', mean(eps_nn(:, pass))*1e9);
    fprintf('  NN residual (detrended) %.2f ns, predicted sigma_eps*sqrt(N/15) = %.2f ns\n', mean(res(:, pass, 2))*1e9, mean(eps_nn(:, pass))*sqrt(N/15)*1e9);
    fprintf('  LS residual, lags 1..%d: %.2f ns (%.1f times smaller than NN)\n', Lmax, mean(res(:, pass, 3))*1e9, mean(res(:, pass, 2))/mean(res(:, pass, 3)));
    fprintf('  Strehl after correction: NN %.3f, LS %.3f\n', mean(S_after(:, pass, 1)), mean(S_after(:, pass, 2)));
end

if isnr == numel(snr_list), drift_noisy = drift/R; end
end
figure('Position', [50 50 900 330]);
subplot(1, 2, 1); plot(1:N, sqrt(drift_noisy)*1e9, 'r'); xlabel('element'); ylabel('rms error (ns)');
title('NN error after detrending, 0 dB SNR, pass 1');
subplot(1, 2, 2); plot(xe*1e3, detr(tau)*1e9, 'k', xe*1e3, detr(tnn)*1e9, 'r', xe*1e3, detr(tls)*1e9, 'b--');
xlabel('element position (mm)'); ylabel('delay (ns)'); title('Last realization, pass 2'); legend('true', 'NN', 'LS');
