%% m05_least_squares_conditioning.m  Module 05: least squares, conditioning, noise amplification
% Purpose : fit an attenuation slope (dB vs MHz) with A\b, pinv and QR; show why
%           centering and scaling the columns matters; show that the normal equations
%           square kappa; tabulate kappa vs column angle; Monte Carlo the noise
%           amplification along the weak singular direction; Tikhonov teaser.
% Change  : sig_dB (spectral noise), the outlier, phi_deg (column angle), lambda grid.
% Expect  : slope ~ -1 dB/MHz, SE = 0.0666 dB/MHz, outlier shift +0.203 dB/MHz;
%           kappa = cot(phi/2): 57.3 at 2 deg, 11.43 at 10 deg; Monte Carlo std along
%           v2 within ~3 percent of sigma/sigma_2; order-7 raw-f fit: relative error ~1e-7
%           with QR, ~1e-2 with the normal equations; Tikhonov at lambda = sigma_2^2
%           halves the v2 noise (0.08 -> 0.04).
% Runtime : ~2 s.  Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);

%% 1. Spectral-difference attenuation model
% log-ratio of power spectra between two depths, dB:  y(f) = y0 - 2*alpha0*dz*f
alpha0 = 0.5;            % dB/cm/MHz (one-way, amplitude)
dz = 1.0;                % cm between the two windows
f = linspace(2, 8, 25).';                % MHz, 25 frequencies
ytrue = 2 - 2*alpha0*dz*f;               % slope -1 dB/MHz (round trip: factor 2)
sig_dB = 0.6;                            % noise std per frequency, dB
y = ytrue + sig_dB*randn(size(f));

A = [ones(size(f)) f];                   % columns: intercept, slope
c_bs   = A \ y;                          % QR-based least squares (what backslash does for tall A)
c_pinv = pinv(A) * y;                    % SVD-based pseudo-inverse
[Q, R] = qr(A, 0);                       % economy QR: A = Q*R, Q'*Q = I
c_qr   = R \ (Q' * y);
c_ne   = (A'*A) \ (A'*y);                % normal equations (fine here, kappa small)
fprintf('slope: \\ %.4f | pinv %.4f | QR %.4f | normal eq %.4f dB/MHz (truth -1)\n', ...
        c_bs(2), c_pinv(2), c_qr(2), c_ne(2));
fprintf('alpha0 estimate = %.3f dB/cm/MHz (truth %.2f)\n', -c_bs(2)/(2*dz), alpha0);
r = y - A*c_bs;
fprintf('residual orthogonal to the columns: max|A''*r| = %.1e\n', max(abs(A'*r)));
SE = sig_dB / sqrt(sum((f - mean(f)).^2));
fprintf('slope standard error sigma/sqrt(sum (f-fbar)^2) = %.4f dB/MHz\n', SE);
Cx = sig_dB^2 * inv(A'*A);               % covariance of the estimate (inv is fine for a 2x2 readout)
fprintf('same from sqrt(sigma^2 (A''A)^-1)(2,2) = %.4f\n', sqrt(Cx(2,2)));

%% 2. One outlier (+6 dB at 7.75 MHz)
yo = y; [~, io] = min(abs(f - 7.75)); yo(io) = yo(io) + 6;
co = A \ yo;
fprintf('outlier at %.2f MHz moves the slope by %+.3f dB/MHz (= %.1f SE); theory 6*(f_i - fbar)/Sxx = %+.3f\n', ...
        f(io), co(2) - c_bs(2), (co(2) - c_bs(2))/SE, 6*(f(io) - mean(f))/sum((f - mean(f)).^2));

%% 3. Column scaling: raw f versus centered t = (f-5)/3
At = [ones(size(f)) (f - 5)/3];
fprintf('cond([1 f]) = %.2f, cond([1 (f-5)/3]) = %.2f\n', cond(A), cond(At));
for p = [3 5 7]
    Araw = zeros(numel(f), p+1); Acen = Araw;
    for k = 0:p, Araw(:,k+1) = f.^k; Acen(:,k+1) = ((f - 5)/3).^k; end
    fprintf('order %d: cond raw = %.1e, centered = %.1f\n', p, cond(Araw), cond(Acen));
end

%% 4. Normal equations square kappa: order-7 polynomial in raw f, noise-free data
p = 7; Araw = zeros(numel(f), p+1);
for k = 0:p, Araw(:,k+1) = f.^k; end
ct = ones(p+1, 1); bt = Araw * ct;       % exact data, known answer
c1 = Araw \ bt;                          % QR
c2 = (Araw'*Araw) \ (Araw'*bt);          % normal equations
fprintf('order 7 raw: cond(A) = %.1e, cond(A''A) = %.1e (= cond(A)^2 = %.1e)\n', cond(Araw), cond(Araw'*Araw), cond(Araw)^2);
fprintf('  relative error: QR %.1e, normal equations %.1e\n', norm(c1 - ct)/norm(ct), norm(c2 - ct)/norm(ct));
fprintf('  digits lost ~ log10(kappa): QR %.1f, normal eq %.1f\n', log10(cond(Araw)), log10(cond(Araw)^2));

%% 5. kappa versus the angle between two unit columns
phis = [90 30 10 5 2];
fprintf('phi(deg)  kappa   cot(phi/2)\n');
for ph = phis
    Ap = [1 cosd(ph); 0 sind(ph)];
    fprintf('%6g  %7.2f  %7.2f\n', ph, cond(Ap), cotd(ph/2));
end

%% 6. Monte Carlo: noise is amplified only along the weak direction v2
phi_deg = 10; sig = 0.01; M = 4000;
Ap = [1 cosd(phi_deg); 0 sind(phi_deg)];
[U, S, V] = svd(Ap); s = diag(S);
xt = [1; 1]; bp = Ap*xt;
Xhat = Ap \ (bp + sig*randn(2, M));     % M solves at once
E = Xhat - xt;
std_v = std(V.' * E, 0, 2);              % std of the error along v1 and v2
fprintf('phi = %g deg: sigma = [%.4f %.4f], std along v1 = %.4f (theory %.4f), v2 = %.4f (theory %.4f)\n', ...
        phi_deg, s, std_v(1), sig/s(1), std_v(2), sig/s(2));

%% 7. Tikhonov teaser: (A'A + lambda I)^-1 A'b
% The bias is lambda/(sigma_i^2 + lambda) times the component of the truth along v_i.
% x = (1,1) lies exactly along v1 (equal unit columns), so shrinking v2 costs it no bias;
% x = (1,0) has a v2 component, so the same lambda buys less.
lams = [0 1e-4 1e-3 s(2)^2 1e-1 1];
truths = [1 1; 1 0].';
rmse = zeros(2, numel(lams));
for q = 1:2
    xq = truths(:,q); bq = Ap*xq; Nq = sig*randn(2, M);
    fprintf('truth x = (%g, %g): v2''x = %+.3f\n', xq, V(:,2).'*xq);
    fprintf('  lambda     f2      std_v2   bias    RMSE\n');
    for k = 1:numel(lams)
        lam = lams(k);
        Mreg = (Ap'*Ap + lam*eye(2)) \ Ap';
        Er = Mreg * (bq + Nq) - xq;
        f2 = s(2)^2/(s(2)^2 + lam);
        rmse(q,k) = sqrt(mean(sum(Er.^2, 1)));
        fprintf('  %8.1e  %5.3f  %7.4f  %6.4f  %6.4f\n', lam, f2, std(V(:,2).'*Er), norm(mean(Er, 2)), rmse(q,k));
    end
end

%% 8. Figures
figure;
subplot(2,2,1);
plot(f, y, 'ko', f, A*c_bs, 'b-', f, ytrue, 'k--', f(io), yo(io), 'rs', f, A*co, 'r-');
xlabel('frequency (MHz)'); ylabel('log-spectral ratio (dB)'); title('LS fit; red: with outlier');
subplot(2,2,2);
plot(E(1,:), E(2,:), '.', 'markersize', 2); axis equal;
hold on; tt = linspace(0, 2*pi, 200);
el = V * diag(2*sig./s) * [cos(tt); sin(tt)]; plot(el(1,:), el(2,:), 'r-', 'linewidth', 1.5); hold off;
xlabel('error in x_1'); ylabel('error in x_2'); title(sprintf('error cloud, phi = %g deg', phi_deg));
subplot(2,2,[3 4]);
semilogx(max(lams, 1e-5), rmse(1,:), 'o-', max(lams, 1e-5), rmse(2,:), 's-');
xlabel('\lambda (0 plotted at 1e-5)'); ylabel('RMSE of x'); legend('x = (1,1)', 'x = (1,0)');
title('Tikhonov RMSE');
