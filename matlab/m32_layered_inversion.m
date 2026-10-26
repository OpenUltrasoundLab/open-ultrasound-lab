%% m32_layered_inversion.m
% From average to local sound speed in a layered medium (straight-ray model).
% Purpose : A global estimator (Section 32.2) returns the average speed down
%           to each depth, c_avg(z) = z / T(z), T(z) = int_0^z s(z') dz'.
%           In slowness this is linear: z * s_avg(z) = int_0^z s dz', so the
%           local slowness is the derivative of z * s_avg. Differencing two
%           noisy averages amplifies their noise by about sqrt(2) z / dz.
%           Compare three inversions of the same noisy c_avg(z):
%             (1) naive differencing on the measurement grid,
%             (2) Tikhonov with a first-difference penalty, lambda chosen by
%                 the discrepancy principle (residual = expected noise),
%             (3) layer-constrained least squares (interfaces known from the
%                 B-mode image, one unknown speed per layer).
% Change  : sig_c (m/s noise on each average), dz (depth step), the profile.
% Expect  : naive local error ~ 1.41 z/dz * sig_c (about 85 m/s at 30 mm for
%           sig_c = 2 m/s, dz = 1 mm); Tikhonov recovers the liver to a few
%           m/s but rounds the interfaces; the layer fit gives each layer to
%           +-1 to 3 m/s, worst for the thin muscle layer. Monte Carlo exceeds
%           the linear noise rule where errors reach ~10 % (1/s is nonlinear).
% Runtime : about 1 s in Octave.
clear; close all; rng(1);

dz = 1e-3; zg = (dz:dz:60e-3)';              % cell bottoms (m), 60 cells of 1 mm
n = numel(zg); zc = zg - dz/2;               % cell centres
iface = [8e-3 14e-3];                        % fat | muscle | liver interfaces (m)
cl = [1450 1585 1550];                       % fat, muscle, liver speeds (m/s)
lay = 1 + (zc > iface(1)) + (zc > iface(2)); % layer index of every cell
c_true = cl(lay)'; s_true = 1./c_true;       % local speed and slowness

km = (4:60)'; zk = zg(km); m = numel(km);    % averages measured at 4..60 mm
A = double(bsxfun(@le, (1:n), km))*dz;       % b_k = z_k * s_avg(z_k) = sum_{j<=k} s_j dz
cavg_true = zk./(A*s_true);
sig_c = 2;                                   % noise on each c_avg (m/s)
cavg = cavg_true + sig_c*randn(m, 1);
b = zk./cavg;                                % measured z * s_avg (s)

% (1) naive differencing
s_naive = diff(b)/dz; z_naive = zk(2:end) - dz/2; c_naive = 1./s_naive;

% predicted noise amplification of (1): sigma_s,k = sigma_savg sqrt(z_k^2 + z_{k-1}^2)/dz
sig_savg = sig_c./cavg_true.^2;
pred = sqrt(zk(2:end).^2 + zk(1:end-1).^2)/dz .* sig_savg(2:end) .* c_true(km(2:end)).^2;   % in m/s

% Monte Carlo check of (1)
R = 400; E = zeros(R, m-1);
for r = 1:R
    bb = zk./(cavg_true + sig_c*randn(m, 1));
    E(r, :) = (1./(diff(bb)/dz))' - c_true(km(2:end))';
end
mc = std(E)';
k30 = find(abs(z_naive - 29.5e-3) < 1e-6); k50 = find(abs(z_naive - 49.5e-3) < 1e-6);
fprintf('Naive differencing, std of local c: predicted %.1f / MC %.1f m/s at 30 mm; %.1f / %.1f at 50 mm (rule 1.41 z/dz sig_c: %.0f, %.0f)\n', ...
    pred(k30), mc(k30), pred(k50), mc(k50), sqrt(2)*30*sig_c, sqrt(2)*50*sig_c);

% (2) Tikhonov with first-difference penalty, discrepancy principle
D = diff(eye(n));
W = diag(1./(zk.*sig_savg));                 % whiten: std of b_k = z_k * sigma_savg
lams = logspace(1, 9, 120); res = zeros(size(lams));
for q = 1:numel(lams)
    s_q = (A'*W'*W*A + lams(q)^2*(D'*D)) \ (A'*W'*W*b);
    res(q) = norm(W*(A*s_q - b))^2;
end
qd = find(res <= m, 1, 'last'); lam = lams(qd);   % largest lambda with chi^2 <= number of data
s_tik = (A'*W'*W*A + lam^2*(D'*D)) \ (A'*W'*W*b); c_tik = 1./s_tik;
liv = zc > 20e-3;
fprintf('Tikhonov (discrepancy, lambda = %.2g): liver (20-60 mm) mean %.1f m/s, RMSE %.1f m/s; whole-profile RMSE %.1f m/s\n', ...
    lam, mean(c_tik(liv)), sqrt(mean((c_tik(liv) - c_true(liv)).^2)), sqrt(mean((c_tik - c_true).^2)));

% (3) layer-constrained least squares: s = P * s_layer, P(j,l) = 1 if cell j is in layer l
P = double(bsxfun(@eq, lay, 1:3));
B = W*A*P; s_lay = B \ (W*b); C = inv(B'*B);  % covariance of the layer slownesses
c_lay = 1./s_lay; sd_lay = sqrt(diag(C)).*c_lay.^2;
fprintf('Layer fit: fat %.1f +- %.1f, muscle %.1f +- %.1f, liver %.1f +- %.1f m/s (true %d, %d, %d)\n', ...
    c_lay(1), sd_lay(1), c_lay(2), sd_lay(2), c_lay(3), sd_lay(3), cl);
fprintf('Average at 60 mm: %.1f m/s, i.e. %.1f m/s below the liver value it is often mistaken for\n', ...
    cavg_true(end), cl(3) - cavg_true(end));

figure('Position', [40 40 1050 430]);
subplot(1, 2, 1);
plot(c_true, zc*1e3, 'k-', 'LineWidth', 2); hold on;
plot(cavg, zk*1e3, 'b.', 'MarkerSize', 10);
plot(c_naive, z_naive*1e3, 'Color', [0.8 0.5 0.5]);
plot(c_tik, zc*1e3, 'r-', 'LineWidth', 1.6);
stairs([c_lay(lay)' c_lay(3)], [0; zg]*1e3, 'g--', 'LineWidth', 1.4);
set(gca, 'YDir', 'reverse'); xlim([1300 1800]); ylim([0 60]); grid on;
xlabel('sound speed (m/s)'); ylabel('depth (mm)');
legend('true local', 'measured average', 'naive difference', 'Tikhonov', 'layer fit', 'Location', 'northeast');
title(sprintf('Average-to-local inversion, sigma_c = %g m/s', sig_c));
subplot(1, 2, 2);
plot(z_naive*1e3, mc, 'k.', z_naive*1e3, pred, 'r-', 'LineWidth', 1.4); grid on;
xlabel('depth (mm)'); ylabel('std of local c (m/s)'); legend('Monte Carlo', 'sqrt(z_k^2 + z_{k-1}^2)/dz rule', 'Location', 'northwest');
title('Noise amplification of naive differencing');
