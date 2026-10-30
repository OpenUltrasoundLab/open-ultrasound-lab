%% m36_tolerance_limit.m  Statistical maximum of acoustic output: one-sided normal tolerance limit
% Purpose : FDA's ultrasound guidance (Feb 2023, section 5.2.4.1) recommends reporting a
%           statistical maximum from the one-sided tolerance limit for normal distributions,
%           with confidence g and population proportion P of 0.9 unless justified:
%               x_stat = xbar + k*s,   k = t'^{-1}_{n-1, delta}(g) / sqrt(n),  delta = z_P*sqrt(n)
%           This script computes k exactly (numerical integration of the noncentral t),
%           compares it with the classical closed-form approximation, checks coverage by
%           Monte Carlo, and applies it to a worked MI / ISPTA.3 example with the
%           guidance's uncertainty rule (+-15 % MI, +-30 % intensities).
% Change  : g, P, nList, the example data (miBar, miS, isBar, isS), U_mi.
% Expected: k(g=P=0.9) = 4.2582 (n=3), 2.7423 (5), 2.0657 (10), 1.7652 (20), 1.6571 (30);
%           MC coverage of xbar+k*s close to 0.90; sample max of 10 covers ~0.65.
% Base MATLAB / GNU Octave only. Runtime a few seconds (about 15 s in Octave).
clear; close all; rng(1);

%% Helpers (anonymous functions only)
Phi  = @(x) 0.5*erfc(-x/sqrt(2));          % standard normal CDF
zq   = @(p) sqrt(2)*erfinv(2*p - 1);       % standard normal quantile
% CDF of the noncentral t with nu dof and noncentrality d, at t:
%   F = int_0^inf Phi(t*w/sqrt(nu) - d) f_chi(w; nu) dw   (w = sqrt(chi-square))
W       = @(nu) linspace(1e-9, sqrt(nu) + 12, 4001);                 % integration grid in w
fchi    = @(w, nu) exp((nu - 1)*log(w) - w.^2/2 - (nu/2 - 1)*log(2) - gammaln(nu/2));
nctcdf0 = @(t, nu, d) trapz(W(nu), Phi(t*W(nu)/sqrt(nu) - d) .* fchi(W(nu), nu));

%% Exact k versus n, and the approximation
g = 0.90; P = 0.90;
nList = [2 3 4 5 6 8 10 15 20 30 50];
zP = zq(P); zg = zq(g);
kEx = zeros(size(nList)); kAp = zeros(size(nList));
for i = 1:numel(nList)
  n = nList(i); nu = n - 1; d = zP*sqrt(n);
  lo = 0; hi = 80;                          % bisection on t for F(t) = g
  for it = 1:60
    mid = 0.5*(lo + hi);
    if nctcdf0(mid, nu, d) < g, lo = mid; else, hi = mid; end
  end
  kEx(i) = 0.5*(lo + hi)/sqrt(n);
  a = 1 - zg^2/(2*(n - 1)); b = zP^2 - zg^2/n;
  kAp(i) = (zP + sqrt(zP^2 - a*b))/a;
end
fprintf('One-sided normal tolerance factor, g = %.2f, P = %.2f (z_P = %.4f)\n', g, P, zP);
fprintf('   n    k exact   k approx   approx error\n');
for i = 1:numel(nList)
  fprintf('%4d   %8.4f   %8.4f   %+6.1f %%\n', nList(i), kEx(i), kAp(i), 100*(kAp(i)/kEx(i) - 1));
end

%% Monte Carlo coverage check (n = 10)
n = 10; k10 = kEx(nList == 10); trials = 20000;
x = randn(n, trials);                      % unit normal production population
xbar = mean(x, 1); s = std(x, 0, 1);
x90 = zP;                                  % true 90th percentile
covTol = mean(xbar + k10*s >= x90);
covMax = mean(max(x, [], 1) >= x90);
fprintf('\nMonte Carlo, n = %d, %d trials:\n', n, trials);
fprintf('  P(xbar + k s >= true P-quantile) = %.3f   (target g = %.2f)\n', covTol, g);
fprintf('  P(sample max  >= true P-quantile) = %.3f   (theory 1 - 0.9^%d = %.3f)\n', covMax, n, 1 - 0.9^n);

%% Worked example: MI and ISPTA.3 statistical maxima, uncertainty rule
miBar = 1.62; miS = 0.06;                  % MI over n probes (illustrative)
isBar = 420;  isS = 45;                    % ISPTA.3 in mW/cm^2 (illustrative)
k3 = kEx(nList == 3);
fprintf('\nWorked example (illustrative data):\n');
fprintf('  MI:      n=10 -> %.3f ; n=3 -> %.3f   (Track 3 limit 1.9, or ISPPA.3 <= 190 W/cm^2)\n', miBar + k10*miS, miBar + k3*miS);
fprintf('  ISPTA.3: n=10 -> %.0f ; n=3 -> %.0f mW/cm^2   (limit 720)\n', isBar + k10*isS, isBar + k3*isS);
U_mi = 20;                                 % your expanded uncertainty on MI, percent
excess = max(0, U_mi - 15);
fprintf('  MI uncertainty %d %% vs FDA estimate 15 %%: limit reduced by %d %% -> %.3f\n', U_mi, excess, 1.9*(1 - excess/100));

%% Figure
figure('Color', 'w');
subplot(1, 2, 1);
plot(nList, kEx, 'o-', 'LineWidth', 1.5); hold on;
plot(nList, kAp, 's--', 'LineWidth', 1.2);
xlabel('units measured n'); ylabel('tolerance factor k');
title('k(n), g = P = 0.9'); legend('exact (noncentral t)', 'approximation'); grid on;
subplot(1, 2, 2);
bar([covTol covMax]); set(gca, 'XTickLabel', {'xbar+ks', 'max of 10'});
hold on; plot([0.5 2.5], [g g], 'k--');
ylabel('coverage of the true 90th percentile'); title('Monte Carlo, n = 10'); ylim([0 1]); grid on;
