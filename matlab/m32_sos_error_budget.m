%% m32_sos_error_budget.m
% Error budget of a global sound-speed estimate from wavefront curvature:
% how per-element delay jitter becomes m/s (Fisher information / CRLB).
% Purpose : A point (or a coherent speckle patch) at depth z returns echoes
%           whose receive arrival times across the aperture follow
%              t_i = t_a + ( sqrt(z^2 + x_i^2) - z ) / c,   z = c t_a,
%           with t_a the apex (on-axis) one-way time, measured from the echo.
%           Unknowns: c and an arrival offset t0 (common delay). The
%           parabolic (small-aperture) limit gives a closed form
%              sigma_c = c^2 z sigma_t sqrt(45/(4N)) / X^2
%                      = 4 c^2 F^2 sigma_t sqrt(45/(4N)) / z,   F = z/(2X),
%           for N elements uniformly spread over [-X, X] and independent
%           per-element jitter sigma_t. (A) checks it against the exact
%           Fisher matrix of the hyperbolic model and a Monte Carlo
%           least-squares fit; (B) tabulates sigma_c against depth and
%           aperture; (C) shows the bias from ignoring a 1 ns/element
%           systematic delay error (e.g. lens or element delay calibration).
% Change  : sigma_t, N, pitch, z_list, the systematic error pattern.
% Expect  : closed form, exact CRLB and Monte Carlo agree within ~5 % (1.62,
%           1.69, 1.68 m/s for 64 elements, 5 ns jitter, z = 30 mm);
%           sigma_c grows linearly with z at a fixed aperture and falls as
%           N^-2.5 at fixed pitch (aperture^2 times sqrt(N)); a 1 ns
%           quadratic calibration error at the aperture edge biases c by
%           about -0.8 m/s (0.05 % of c per ns of the ~970 ns curvature),
%           comparable to the CRLB and, unlike jitter, it never averages out.
% Runtime : about 2 s in Octave.
clear; close all; rng(3);

c = 1540; pitch = 0.3e-3;                    % true speed (m/s), pitch (m)
sigma_t = 5e-9;                              % per-element arrival-time jitter (s)
model = @(cc, t0, x, ta) t0 + ta + (sqrt((cc*ta)^2 + x.^2) - cc*ta)/cc;   % one-way receive times (s)

% (A) one configuration: 64 elements, z = 30 mm
N = 64; x = ((0:N-1) - (N-1)/2)*pitch; X = max(abs(x)); z = 30e-3; ta = z/c;
h = 1e-3;                                     % finite-difference step for the Jacobian (m/s and s)
J = [(model(c + h, 0, x, ta) - model(c - h, 0, x, ta))'/(2*h), ones(N, 1)];   % d t / d(c, t0)
F = (J'*J)/sigma_t^2; Ci = inv(F); crlb = sqrt(Ci(1, 1));
Xrms = sqrt(mean(x.^4) - mean(x.^2)^2);      % discrete spread of x^2 (continuous: X^2 sqrt(4/45))
closed = c^2*z*sigma_t/(sqrt(N)*Xrms);
closed_u = c^2*z*sigma_t*sqrt(45/(4*N))/X^2;
R = 3000; est = zeros(R, 1);
for r = 1:R
    tm = model(c, 0, x, ta) + sigma_t*randn(1, N);
    p = [c; 0];                              % Gauss-Newton on (c, t0)
    for it = 1:6
        res = tm' - model(p(1), p(2), x, ta)';
        Jg = [(model(p(1) + h, p(2), x, ta) - model(p(1) - h, p(2), x, ta))'/(2*h), ones(N, 1)];
        p = p + Jg \ res;
    end
    est(r) = p(1);
end
fprintf('(A) N = %d, aperture %.1f mm, z = %.0f mm (F = %.2f), sigma_t = %.0f ns\n', N, 2*X*1e3, z*1e3, z/(2*X), sigma_t*1e9);
fprintf('    exact CRLB %.2f m/s | parabolic closed form %.2f (uniform-aperture form %.2f) | Monte Carlo std %.2f, bias %+.2f m/s\n', ...
    crlb, closed, closed_u, std(est), mean(est) - c);

% (B) table: sigma_c versus depth and aperture (exact CRLB)
z_list = [10 20 30 40 60 80]*1e-3; N_list = [32 64 128];
fprintf('\n(B) CRLB of c (m/s), sigma_t = %.0f ns\n   z (mm) ', sigma_t*1e9); fprintf('   N = %3d (F)  ', N_list); fprintf('\n');
S = zeros(numel(z_list), numel(N_list));
for iz = 1:numel(z_list)
    fprintf('   %5.0f  ', z_list(iz)*1e3);
    for in = 1:numel(N_list)
        Nn = N_list(in); xx = ((0:Nn-1) - (Nn-1)/2)*pitch; tz = z_list(iz)/c;
        Jn = [(model(c + h, 0, xx, tz) - model(c - h, 0, xx, tz))'/(2*h), ones(Nn, 1)];
        Ci = inv(Jn'*Jn)*sigma_t^2; S(iz, in) = sqrt(Ci(1, 1));
        fprintf('   %7.2f (%4.2f)', S(iz, in), z_list(iz)/(2*max(abs(xx))));
    end
    fprintf('\n');
end

% (C) a systematic error: quadratic delay error of 1 ns at the aperture edge
e_sys = 1e-9*(x/X).^2;                       % s, e.g. lens/element delay calibration
p = [c; 0]; tm = model(c, 0, x, ta) + e_sys;
for it = 1:8
    res = tm' - model(p(1), p(2), x, ta)';
    Jg = [(model(p(1) + h, p(2), x, ta) - model(p(1) - h, p(2), x, ta))'/(2*h), ones(N, 1)];
    p = p + Jg \ res;
end
fprintf('\n(C) 1 ns quadratic calibration error at the aperture edge: bias %+.2f m/s (CRLB at 5 ns jitter: %.2f m/s)\n', p(1) - c, crlb);
fprintf('    rule: the curvature of a 1540 m/s wavefront across this aperture is X^2/(2 z c) = %.0f ns; 1 ns of it is %.2f %% of c\n', ...
    X^2/(2*z*c)*1e9, 100*1e-9/(X^2/(2*z*c))/2);

figure('Position', [60 60 980 380]);
subplot(1, 2, 1); hist(est, 40); xlabel('estimated c (m/s)'); ylabel('count');
title(sprintf('Monte Carlo, std %.2f m/s, CRLB %.2f m/s', std(est), crlb)); grid on;
subplot(1, 2, 2); loglog(z_list*1e3, S, 'o-', 'LineWidth', 1.4); grid on;
xlabel('depth z (mm)'); ylabel('CRLB of c (m/s)'); legend('N = 32', 'N = 64', 'N = 128', 'Location', 'northwest');
title(sprintf('Wavefront-curvature CRLB, sigma_t = %.0f ns', sigma_t*1e9));
