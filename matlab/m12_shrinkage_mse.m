%% m12_shrinkage_mse.m
% Purpose : biased estimators can have lower MSE than the CRLB, and the
%           Bayesian view explains when that is legitimate.
%   Part 1: scalar shrinkage theta_hat = c * ybar, ybar ~ N(theta, s^2),
%           s^2 = sigma^2/N (the CRLB). MSE = c^2 s^2 + (1-c)^2 theta^2,
%           below the CRLB only for |theta| < s*sqrt((1+c)/(1-c)).
%   Part 2: James-Stein for p >= 3 means, y ~ N(theta, I_p):
%           theta_JS = (1 - (p-2)/||y||^2) y. Exact risk
%           R = p - (p-2)^2 E[1/(p-2+2K)], with K the Poisson mixing index
%           (mean ||theta||^2/2) of the noncentral chi-square ||y||^2,
%           always < p (the MLE risk). Positive-part version is better still.
%   Part 3: Bayes. Prior theta ~ N(0, tau^2): the posterior mean is
%           shrinkage with c = tau^2/(tau^2 + s^2); its Bayesian MSE equals
%           the Bayesian CRB 1/(1/s^2 + 1/tau^2) (data + prior information).
% Change  : c, p, tau, the theta grid.
% Expect  : Part 1, c = 0.8: break-even |theta| = 3.00 s. Part 2: risk 2.00
%           at theta = 0 for every p; MC matches the exact formula to ~1 %.
%           Part 3: Bayesian MSE = Bayesian CRB (ratio 1.00).
rng(1);

%% Part 1: scalar shrinkage
s = 1;                                   % CRLB std of the unbiased estimate
cList = [0.5 0.8 0.95];
th = linspace(-5, 5, 401);
fprintf('Part 1: scalar shrinkage c*ybar, MSE in units of the CRLB\n');
for c = cList
    fprintf('  c = %.2f: MSE(0) = %.3f, break-even |theta| = %.2f s, MSE(3 s) = %.2f\n', ...
        c, c^2, sqrt((1 + c)/(1 - c)), c^2 + (1 - c)^2*9);
end
M = 200000; theta0 = 1.0; c = 0.8;
yb = theta0 + s*randn(M, 1);
fprintf('  Monte Carlo, theta = %.1f s, c = %.1f: MSE = %.4f (theory %.4f)\n', theta0, c, ...
    mean((c*yb - theta0).^2), c^2 + (1 - c)^2*theta0^2);

%% Part 2: James-Stein
jsRisk = @(p, lam) p - (p - 2)^2 * sum(exp(-lam/2 + (0:400)*log(max(lam/2, realmin)) - gammaln(1:401)) ./ (p - 2 + 2*(0:400)));
pList = [3 10 50]; normList = [0 1 2 3];   % ||theta|| / sqrt(p)
fprintf('\nPart 2: James-Stein risk (MLE risk = p)\n');
fprintf('     p  ||th||/sqrt(p)   exact     MC JS   MC JS+\n');
riskCurve = zeros(numel(pList), 31); nn = linspace(0, 3, 31);
for ip = 1:numel(pList)
    p = pList(ip);
    for k = 1:numel(normList)
        theta = zeros(1, p); theta(1) = normList(k)*sqrt(p);
        Y = theta + randn(40000, p);           % implicit expansion (R2016b+, Octave)
        n2 = sum(Y.^2, 2);
        JS = (1 - (p - 2)./n2) .* Y; JSp = max(0, 1 - (p - 2)./n2) .* Y;
        fprintf('  %4d  %10.1f  %9.3f %9.3f %8.3f\n', p, normList(k), jsRisk(p, sum(theta.^2)), ...
            mean(sum((JS - theta).^2, 2)), mean(sum((JSp - theta).^2, 2)));
    end
    for k = 1:numel(nn), riskCurve(ip, k) = jsRisk(p, nn(k)^2*p)/p; end
end

%% Part 3: Bayes posterior mean = shrinkage; Bayesian CRB
tau = 1.5; cB = tau^2/(tau^2 + s^2);
thetaB = tau*randn(M, 1); yB = thetaB + s*randn(M, 1);
mseB = mean((cB*yB - thetaB).^2); bcrb = 1/(1/s^2 + 1/tau^2);
fprintf('\nPart 3: prior std tau = %.1f s: c = %.3f, Bayesian MSE = %.4f, Bayesian CRB = %.4f, ratio %.3f\n', ...
    tau, cB, mseB, bcrb, mseB/bcrb);
fprintf('  same estimator, MSE at theta = 3 tau = %.2f (CRLB 1): the prior must be right\n', ...
    cB^2 + (1 - cB)^2*(3*tau)^2);

%% Figures
figure('Position', [80 80 950 360]);
subplot(1, 2, 1);
for c = cList, plot(th, c^2 + (1 - c)^2*th.^2, 'LineWidth', 1.3); hold on; end
plot(th, ones(size(th)), '--', 'Color', [0.6 0 0], 'LineWidth', 2); ylim([0 3]); grid on;
xlabel('\theta / s  (s = CRLB std)'); ylabel('MSE / CRLB');
title('Shrinkage c\cdot y: below the CRLB near \theta = 0 only');
legend('c = 0.5', 'c = 0.8', 'c = 0.95', 'unbiased (CRLB)', 'Location', 'north');
subplot(1, 2, 2);
plot(nn, riskCurve.', 'LineWidth', 1.3); hold on; plot(nn, ones(size(nn)), '--', 'Color', [0.6 0 0], 'LineWidth', 2);
xlabel('||\theta|| / sqrt(p)'); ylabel('risk / p'); ylim([0 1.1]); grid on;
title('James-Stein dominates the MLE for p >= 3');
legend('p = 3', 'p = 10', 'p = 50', 'MLE', 'Location', 'southeast');
