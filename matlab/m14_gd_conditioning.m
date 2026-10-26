%% m14_gd_conditioning.m  --  gradient descent, momentum and the condition number
% Purpose : minimize the quadratic f(x) = 1/2 x'*A*x (minimum at 0) with eigenvalues
%           mu = 1 and L = kappa, and count iterations to reduce the error norm by 1e-6.
%           Compares gradient descent with the optimal fixed step 2/(L+mu), with step 1/L,
%           Polyak heavy-ball with its optimal (alpha, beta), and Nesterov's method
%           (step 1/L, momentum (sqrt(kappa)-1)/(sqrt(kappa)+1)) against the rate formulas.
%           Then shows that diagonal rescaling of badly scaled parameters (here a sound
%           speed in m/s and a dimensionless amplitude) removes most of kappa.
% Change  : kappa_list, tol, the rotation angle th, the momentum settings.
% Expect  : kappa = 100: GD(2/(L+mu)) 691 iterations (theory 691), GD(1/L) 1341 (1375; this
%           start puts only 1/sqrt(2) of the error on the slow mode), heavy ball 93 and
%           Nesterov 155 (rate-only predictions 69 and 131: both tuned iterations have a
%           double root, so the error decays as k*rho^k, not rho^k). Jacobi scaling of a
%           (m/s, amplitude) problem drops kappa(J'J) from about 1e6 to about 1.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);

%% 1. Iterations to tolerance vs kappa
tol = 1e-6;  kappa_list = [10 100 1000];  maxit = 1e5;
th = pi/6;  Q = [cos(th) -sin(th); sin(th) cos(th)];          % rotate so the axes are not the eigenvectors
fprintf('kappa | GD 2/(L+mu)  theory | GD 1/L  theory | heavy ball  theory | Nesterov  theory\n');
hist100 = {};
for kappa = kappa_list
    mu = 1;  L = kappa;  A = Q * diag([mu L]) * Q.';
    x0 = Q * [1; 1];
    sk = sqrt(kappa);
    cfg = { 2/(L+mu), 0, 'hb';                                  % {alpha, beta, type}
            1/L,      0, 'hb';
            4/(sqrt(L)+sqrt(mu))^2, ((sk-1)/(sk+1))^2, 'hb';
            1/L,      (sk-1)/(sk+1), 'nest' };
    its = zeros(1, 4);
    for m = 1:4
        a = cfg{m,1};  b = cfg{m,2};  x = x0;  xp = x0;  e = zeros(1, maxit);
        for k = 1:maxit
            e(k) = norm(x) / norm(x0);
            if e(k) < tol, break; end
            if strcmp(cfg{m,3}, 'nest')
                y = x + b*(x - xp);  xn = y - a*(A*y);
            else
                xn = x - a*(A*x) + b*(x - xp);
            end
            xp = x;  x = xn;
        end
        its(m) = k - 1;
        if kappa == 100, hist100{m} = e(1:k); end
    end
    th_gd = log(1/tol) / log((kappa+1)/(kappa-1));
    th_1L = log(1/tol) / -log(1 - 1/kappa);
    th_hb = log(1/tol) / log((sk+1)/(sk-1));
    th_ne = log(1/tol) / -log(1 - 1/sk);
    fprintf('%5d | %6d  %6.0f     | %6d  %6.0f | %6d  %6.0f     | %6d  %6.0f\n', ...
        kappa, its(1), th_gd, its(2), th_1L, its(3), th_hb, its(4), th_ne);
end

%% 2. Divergence threshold: step just above 2/L
kappa = 100;  A = Q * diag([1 kappa]) * Q.';  x = Q * [1; 1];
for k = 1:50, x = x - (2.05/kappa) * (A*x); end
fprintf('step 2.05/L after 50 iterations: |x| = %.3g (diverges; factor |1 - 2.05| = 1.05 per step)\n', norm(x));

%% 3. Badly scaled parameters: sound speed (m/s) and amplitude (unitless)
% Linearized model y = J*[dc; dA] with column sensitivities that differ by 3 orders.
J = [randn(50,1) * 1e-3, randn(50,1)];          % d y / d c (per m/s) and d y / d A
H = J.' * J;  D = diag(1 ./ sqrt(diag(H)));     % Jacobi (diagonal) preconditioner
fprintf('kappa(J''J) raw %.3g  ->  after column scaling %.3g\n', cond(H), cond(D*H*D));
kr = cond(H);  ks = cond(D*H*D);
fprintf('GD iterations to 1e-6 (optimal step): raw %.3g, scaled %.0f\n', ...
    log(1e6)/log((kr+1)/(kr-1)), log(1e6)/log((ks+1)/(ks-1)));

%% 4. Figures: error histories (kappa = 100) and paths on the contours
figure('Position', [100 100 1050 420]);
subplot(1,2,1);
cl = {[0 0.3 0.8], [0.6 0.6 0.6], [0.8 0.2 0.5], [0 0.55 0.4]};  hh = zeros(1,6);
for m = 1:4
    hh(m) = semilogy(0:numel(hist100{m})-1, hist100{m}, 'Color', cl{m}, 'LineWidth', 1.5); hold on;
end
kk = 0:700;
hh(5) = semilogy(kk, ((100-1)/(100+1)).^kk, '--', 'Color', cl{1});
hh(6) = semilogy(kk, ((10-1)/(10+1)).^kk, '--', 'Color', cl{3});
grid on; ylim([1e-7 2]); xlim([0 700]); xlabel('iteration k'); ylabel('||x_k|| / ||x_0||');
legend(hh, {'GD 2/(L+\mu)', 'GD 1/L', 'heavy ball', 'Nesterov', '((\kappa-1)/(\kappa+1))^k', '((\surd\kappa-1)/(\surd\kappa+1))^k'}, 'Location', 'northeast');
title('\kappa = 100: measured error vs rate formulas');
subplot(1,2,2);
kappa = 30;  A = Q * diag([1 kappa]) * Q.';
[g1, g2] = meshgrid(linspace(-1.6, 1.6, 200), linspace(-1.6, 1.6, 200));
F = 0.5 * (A(1,1)*g1.^2 + 2*A(1,2)*g1.*g2 + A(2,2)*g2.^2);
contour(g1, g2, log10(F + 1e-3), 14); hold on; colormap(gray);
x0 = Q * [1.2; 0.25];  hp = zeros(1,2);
for m = 1:2
    if m == 1, a = 2/(kappa+1); b = 0; col = [0 0.3 0.8]; else, sk = sqrt(kappa); a = 4/(sk+1)^2; b = ((sk-1)/(sk+1))^2; col = [0.8 0.2 0.5]; end
    x = x0; xp = x0; P = x0;
    for k = 1:60, xn = x - a*(A*x) + b*(x - xp); xp = x; x = xn; P(:, end+1) = x; end
    hp(m) = plot(P(1,:), P(2,:), '.-', 'Color', col, 'LineWidth', 1.2, 'MarkerSize', 10);
end
axis image; xlabel('x_1'); ylabel('x_2'); title('\kappa = 30: GD zig-zags, heavy ball spirals in');
legend(hp, {'GD 2/(L+\mu)', 'heavy ball (optimal \alpha, \beta)'}, 'Location', 'southoutside');
