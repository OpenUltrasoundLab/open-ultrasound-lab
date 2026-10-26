%% m13_tv_vs_tikhonov.m
% Average-to-local sound speed: naive inversion, Tikhonov (L = I, D1, D2) and TV.
% Purpose : A depth-dependent average sound speed cbar(z) (for example from
%           autofocus at each depth, vertical-path approximation) relates to
%           the local slowness s(z) = 1/c(z) by the running mean
%              1/cbar(z) = (1/z) * integral_0^z s(z') dz'      (linear in s).
%           Naive inversion differentiates z/cbar(z), amplifying noise by
%           about sqrt(2) z/dz. Regularize the slowness deviation
%           x = s - 1/c0 (ns/mm) with
%              Tikhonov : min ||A x - b||^2 + lam^2 ||L x||^2,  L = I, D1 or D2
%              TV       : min 0.5||A x - b||^2 + mu ||D1 x||_1  (ADMM, rho = 0.1)
%           Also checks the 1D TV-denoising rule: an isolated plateau of width w
%           (samples) and height 1 comes back with height 1 - 2 mu / w.
% Change  : profile ('layers' or 'ramp'), noise std sc (m/s), dz.
% Expect  : naive std grows linearly with depth (about 110 m/s at 30 mm for
%           1 m/s noise on cbar); TV beats every Tikhonov on the layers and keeps
%           the 8 mm step sharp; on the ramp segment (z > 24 mm) TV staircases,
%           so compare the 'z > 24 mm' column, not only the total RMSE.
% Runtime : about 10 s (the ADMM parameter scan dominates).
clear; close all; rng(1);

dz = 0.4; n = 100; z = dz*(1:n)'; c0 = 1540;      % mm, samples, mm, m/s
A  = tril(ones(n))*dz./repmat(z, 1, n);           % running-mean operator
D1 = diff(eye(n)); D2 = diff(eye(n), 2);
sv = svd(A); fprintf('running-mean operator: sigma_max %.2f, sigma_min %.4f, kappa %.0f\n', sv(1), sv(end), sv(1)/sv(end));
toS = @(c) (1./c - 1/c0)*1e6;                     % m/s -> slowness deviation (ns/mm)
toC = @(x) 1./(1/c0 + x*1e-6);                    % back to m/s
prof = @(kind) 1450*(z <= 8) + 1580*(z > 8 & z <= 20) + 1560*(z > 20 & z <= 24) + ...
    (1560 + strcmp(kind, 'ramp')*2.0*(z - 24)).*(z > 24);

%% Naive inversion: noise growth with depth
c = prof('layers'); cbar = 1./(A*(1./c)); sc = 1.0; M = 2000;
err = zeros(M, n);
for m = 1:M
    sb = 1./(cbar + sc*randn(n, 1));
    sn = [sb(1); (z(2:end).*sb(2:end) - z(1:end-1).*sb(1:end-1))/dz];
    err(m, :) = (1./sn - c)';
end
iz = [25 50 75];
% first-order theory: std(c) = (c/cbar)^2 * sqrt(z_i^2 + z_(i-1)^2)/dz * sc  ~  sqrt(2) z/dz * sc
fprintf('naive inverse, cbar noise %.1f m/s: std of c at z = 10, 20, 30 mm: %s m/s (theory %s)\n', ...
    sc, mat2str(round(std(err(:, iz)))), mat2str(round((c(iz)'./cbar(iz)').^2.*sqrt(z(iz).^2 + z(iz-1).^2)'/dz*sc)));

%% Regularized inversions
edge = @(ce) z(find(ce > 1450 + 0.9*130 & z > 3, 1)) - z(find(ce > 1450 + 0.1*130 & z > 3, 1));  % 10-90 % width of the 8 mm step
hf = figure('Position', [50 50 1000 380]);
kinds = {'layers', 'ramp'}; lamg = logspace(-4, 1, 60); mug = logspace(-3, 1, 25); rho = 0.1;
for kk = 1:2
    c = prof(kinds{kk}); sc = 2.0;
    b = toS(1./(A*(1./c)) + sc*randn(n, 1));
    Ls = {eye(n), D1, D2}; Ln = {'I', 'D1', 'D2'};
    fprintf('\n%s, cbar noise %.1f m/s (oracle parameters, RMSE in m/s):\n', kinds{kk}, sc);
    for q = 1:3
        L = Ls{q}; e = zeros(size(lamg));
        for j = 1:numel(lamg), e(j) = sqrt(mean((toC((A'*A + lamg(j)^2*(L'*L))\(A'*b)) - c).^2)); end
        [eb, jb] = min(e); ce = toC((A'*A + lamg(jb)^2*(L'*L))\(A'*b));
        fprintf('  Tikhonov L = %-2s  lambda %.3g  RMSE %5.1f  (z > 24 mm: %5.1f)  edge 10-90 %% %.1f mm\n', Ln{q}, lamg(jb), eb, sqrt(mean((ce(z > 24) - c(z > 24)).^2)), edge(ce));
        if q == 2, cT = ce; end
    end
    Mi = inv(A'*A + rho*(D1'*D1)); Atb = A'*b; e = zeros(size(mug));
    for j = 1:numel(mug)
        x = Mi*Atb; zz = D1*x; u = zeros(n-1, 1); th = mug(j)/rho;
        for it = 1:1500
            x = Mi*(Atb + rho*D1'*(zz - u)); w = D1*x + u; zz = sign(w).*max(abs(w) - th, 0); u = w - zz;
        end
        e(j) = sqrt(mean((toC(x) - c).^2));
        if j == 1 || e(j) < min(e(1:j-1)), cV = toC(x); end
    end
    [eb, jb] = min(e);
    fprintf('  TV (ADMM)       mu     %.3g  RMSE %5.1f  (z > 24 mm: %5.1f)  edge 10-90 %% %.1f mm\n', mug(jb), eb, sqrt(mean((cV(z > 24) - c(z > 24)).^2)), edge(cV));
    subplot(1, 2, kk);
    plot(z, c, 'k', z, cT, 'm', z, cV, 'b', 'LineWidth', 1.2);
    xlabel('depth z (mm)'); ylabel('local sound speed (m/s)'); title([kinds{kk} ': truth, Tikhonov-D1 (magenta), TV (blue)']);
    ylim([1420 1640]);
end

%% TV denoising contrast-loss rule (A = I)
m = 128; D = diff(eye(m)); y = zeros(m, 1); y(41:60) = 1;   % plateau width w = 20 samples
Mi = inv(eye(m) + 1*(D'*D));
fprintf('\nTV denoising of a 20-sample plateau of height 1:\n');
for mu = [0.5 1 2 4]
    x = y; zz = D*x; u = zeros(m-1, 1);
    for it = 1:3000
        x = Mi*(y + D'*(zz - u)); w = D*x + u; zz = sign(w).*max(abs(w) - mu, 0); u = w - zz;
    end
    fprintf('  mu = %.1f: plateau %.3f (rule 1 - 2mu/w = %.3f), background next to it %.4f\n', mu, x(50), 1 - 2*mu/20, x(10));
end
