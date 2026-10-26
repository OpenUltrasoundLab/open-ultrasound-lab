%% m07_covariance_basics.m  --  covariance ellipse, principal axes, regression lines, Fisher z
% Purpose : draw a 2D Gaussian with a chosen covariance Sigma, estimate the sample
%           covariance, read it as an ellipse (eigenvectors = principal axes), compare
%           the three lines through the cloud, then Monte-Carlo the sampling
%           distribution of the correlation coefficient and check Fisher's z.
% Change  : Sigma (any symmetric positive definite 2x2), Npts, rho_mc, Npairs.
% Expect  : rho = 0.40, lambda = 9.966 / 3.034, axis angle 68.1 deg,
%           slopes y-on-x 0.600, principal 2.488, x-on-y 3.750,
%           39.3 % inside the 1-sigma ellipse, 86.5 % inside 2-sigma;
%           Fisher z std ~ 1/sqrt(N-3) and ~95 % interval coverage.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);

%% 1. Draw samples with covariance Sigma
Sigma = [4 2.4; 2.4 9];          % covariance (a.u.^2), same as Self-check Q1
Npts  = 2000;                    % number of samples
Lc    = chol(Sigma, 'lower');    % Sigma = Lc*Lc'
X     = Lc * randn(2, Npts);     % 2 x Npts, zero mean

%% 2. Sample covariance (1/(N-1)) written out, then rho
mu    = mean(X, 2);
Xc    = X - mu;                              % implicit expansion
S     = (Xc * Xc') / (Npts - 1);             % sample covariance
rho_true = Sigma(1,2) / sqrt(Sigma(1,1) * Sigma(2,2));
rho_hat  = S(1,2) / sqrt(S(1,1) * S(2,2));
fprintf('true Sigma   = [%6.3f %6.3f; %6.3f %6.3f]\n', Sigma');
fprintf('sample S     = [%6.3f %6.3f; %6.3f %6.3f]\n', S');
fprintf('rho true %.3f   rho_hat %.3f\n', rho_true, rho_hat);

%% 3. Eigen-decomposition = principal axes of the ellipse
[U, Lam] = eig(Sigma);
[lam, order] = sort(diag(Lam), 'descend');  U = U(:, order);
if U(1,1) < 0, U(:,1) = -U(:,1); end        % fix sign for a readable angle
theta = atan2(U(2,1), U(1,1)) * 180 / pi;   % major-axis angle from x axis (deg)
fprintf('eigenvalues  = %.3f, %.3f   major axis angle = %.1f deg\n', lam(1), lam(2), theta);

%% 4. The three lines through the cloud
slope_yx  = Sigma(1,2) / Sigma(1,1);        % E[y|x]: minimises vertical errors
slope_pca = U(2,1) / U(1,1);                % principal axis: perpendicular errors (TLS)
slope_xy  = Sigma(2,2) / Sigma(1,2);        % E[x|y] drawn in the x-y plane
fprintf('slopes: y-on-x %.3f   principal %.3f   x-on-y %.3f\n', slope_yx, slope_pca, slope_xy);

%% 5. Fraction inside the k-sigma ellipses (chi-square with 2 dof)
d2 = sum(X .* (Sigma \ X), 1);              % squared Mahalanobis distance
fprintf('inside 1-sigma: %.1f %% (theory %.1f %%)   2-sigma: %.1f %% (theory %.1f %%)\n', ...
    100*mean(d2 <= 1), 100*(1-exp(-1/2)), 100*mean(d2 <= 4), 100*(1-exp(-2)));

figure('Position', [100 100 1000 420]);
subplot(1,2,1);
plot(X(1,1:800), X(2,1:800), '.', 'Color', [0.45 0.6 0.85]); hold on;
t = linspace(0, 2*pi, 200);
for k = [1 2]
    E = U * diag(k*sqrt(lam)) * [cos(t); sin(t)];
    plot(E(1,:), E(2,:), 'k-', 'LineWidth', 1.5);
end
xr = 8;
plot([-xr xr], slope_yx*[-xr xr], 'r-', 'LineWidth', 1.5);
plot([-xr xr], slope_pca*[-xr xr], 'k--', 'LineWidth', 1.2);
plot([-xr xr], slope_xy*[-xr xr], 'Color', [0 0.6 0.5], 'LineWidth', 1.5);
axis equal; xlim([-8 8]); ylim([-10 10]); grid on;
xlabel('x (a.u.)'); ylabel('y (a.u.)');
title('1\sigma/2\sigma ellipses, principal axis (--), y on x (red), x on y (green)');

%% 6. Fisher z Monte Carlo for the sample correlation coefficient
rho_mc = 0.5;  Npairs = 20;  Nexp = 5000;    % true rho, pairs per experiment, experiments
x = randn(Npairs, Nexp);
y = rho_mc * x + sqrt(1 - rho_mc^2) * randn(Npairs, Nexp);
xc = x - mean(x, 1);  yc = y - mean(y, 1);
r  = sum(xc .* yc, 1) ./ sqrt(sum(xc.^2, 1) .* sum(yc.^2, 1));
z  = atanh(r);
se = 1 / sqrt(Npairs - 3);
covered = abs(z - atanh(rho_mc)) <= 1.96 * se;
fprintf('rho = %.2f, N = %d: mean rho_hat %.4f, std %.4f (delta method (1-rho^2)/sqrt(N-1) = %.4f)\n', ...
    rho_mc, Npairs, mean(r), std(r), (1 - rho_mc^2) / sqrt(Npairs - 1));
fprintf('Fisher z: std %.4f (1/sqrt(N-3) = %.4f), 95%% interval coverage %.1f %%\n', std(z), se, 100*mean(covered));

subplot(1,2,2);
edges = linspace(-0.3, 1, 66);  centers = (edges(1:end-1) + edges(2:end)) / 2;  bw = edges(2) - edges(1);
cnt = hist(r, centers);
bar(centers, cnt / (Nexp * bw), 1, 'FaceColor', [0.6 0.7 0.9], 'EdgeColor', 'none'); hold on;
rr = linspace(-0.29, 0.99, 400);
mz = atanh(rho_mc) + rho_mc / (2 * (Npairs - 1));
pdf_r = exp(-0.5 * ((atanh(rr) - mz) / se).^2) ./ (se * sqrt(2*pi) * (1 - rr.^2));
plot(rr, pdf_r, 'r-', 'LineWidth', 1.8);
plot([rho_mc rho_mc], [0 max(pdf_r)*1.1], 'k--');
xlabel('sample correlation coefficient'); ylabel('density');
title(sprintf('\\rho = %.1f, N = %d: histogram vs Fisher-z density', rho_mc, Npairs));
