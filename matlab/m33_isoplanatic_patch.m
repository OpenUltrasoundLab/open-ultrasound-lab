%% m33_isoplanatic_patch.m
% Isoplanatic patch of a phase screen away from the array: how far from the
% guide point one aberration law stays useful, by Monte Carlo against the
% exact Gaussian-screen formula and its large-aperture limit.
% Purpose : 128 elements (pitch 0.3 mm), 5 MHz, focal depth z = 40 mm. A thin
%           Gaussian delay screen (rms sig_tau, correlation rho(dx) =
%           exp(-dx^2/l^2)) lies at depth d. The straight ray from element x_u
%           to the point (x_p, z) crosses it at a_u + x_p d/z, a_u = x_u(1 - d/z).
%           The law measured at x_p = 0 is applied at x_p = dx: the residual is
%           r_u = phi(a_u + delta) - phi(a_u), delta = dx d / z.
%           Exact mean Strehl (Gaussian statistics):
%             E[S] = (1/N^2) sum_uv exp(-sig^2 [2 - 2 rho(delta) - 2 rho(da)
%                                        + rho(da + delta) + rho(da - delta)])
%           Large aperture: S ~ exp(-2 sig^2 (1 - rho(delta))).
%           Correction beats no correction while rho(delta) > 1/2, i.e.
%             |dx| < dx_x = l z sqrt(ln 2) / d.
% Change  : sig_tau, ell, d_list, z.
% Expect  : Monte Carlo on the exact curve; the patch shrinks as z/d; beyond
%           dx_x the borrowed law is worse than none.
% Runtime : about 25 s in Octave.
clear; close all; rng(4);

c = 1540; f0 = 5e6; N = 128; p = 0.3e-3; xe = ((0:N-1) - (N-1)/2)*p; z = 40e-3;
sig_tau = 40e-9; ell = 4e-3; sp = 2*pi*f0*sig_tau; d_list = [5 10 20]*1e-3; R = 200;
rho = @(x) exp(-(x/ell).^2);
dxs = 0.05e-3; xg = -40e-3:dxs:40e-3; g = exp(-2*(((-3*ell):dxs:(3*ell))/ell).^2); g = g/sqrt(sum(g.^2));
dx = [0:0.5:12, 13:1:30]*1e-3;                                                           % lateral distance from the guide point
fprintf('sigma_tau = %.0f ns (%.2f rad), l = %.0f mm, z = %.0f mm; uncorrected E[S] = %.3f\n', sig_tau*1e9, sp, ell*1e3, z*1e3, ...
    sum(sum(exp(-sp^2*(1 - rho(xe' - xe)))))/N^2);
fprintf('\n  d     dx_x = l z sqrt(ln2)/d (exact crossing)   dx(S=1/2) exact   MC S at 0, 2, 4, 8 mm         exact S at 0, 2, 4, 8 mm\n');
Smc = zeros(numel(d_list), numel(dx)); Sex = Smc;
for a = 1:numel(d_list)
    d = d_list(a); au = xe*(1 - d/z); DA = au' - au;
    for q = 1:numel(dx)
        del = dx(q)*d/z;
        Sex(a, q) = sum(sum(exp(-sp^2*(2 - 2*rho(del) - 2*rho(DA) + rho(DA + del) + rho(DA - del)))))/N^2;
    end
    acc = zeros(1, numel(dx));
    for r = 1:R
        phi = 2*pi*f0*sig_tau*conv(randn(1, numel(xg)), g, 'same');
        p0 = interp1(xg, phi, au);
        for q = 1:numel(dx)
            acc(q) = acc(q) + abs(mean(exp(1j*(interp1(xg, phi, au + dx(q)*d/z) - p0))))^2;
        end
    end
    Smc(a, :) = acc/R;
    S0 = sum(sum(exp(-sp^2*(1 - rho(DA)))))/N^2;                              % uncorrected at this d (screen compressed by 1 - d/z)
    half = interp1(Sex(a, :), dx, 0.5); ii = find(Sex(a, :) < S0, 1); cross = interp1(Sex(a, ii-1:ii), dx(ii-1:ii), S0);
    k = arrayfun(@(v) find(abs(dx - v) < 1e-9), [0 2 4 8]*1e-3);
    fprintf('%3.0f mm    %6.2f mm (%5.2f mm)                %6.2f mm         %s   %s  (uncorrected %.3f)\n', d*1e3, ell*z*sqrt(log(2))/d*1e3, cross*1e3, half*1e3, ...
        sprintf('%.3f ', Smc(a, k)), sprintf('%.3f ', Sex(a, k)), S0);
end
fprintf('\nLarge-aperture limit exp(-2 sig^2 (1 - rho)) at delta = l: %.3f; correction crosses the uncorrected level at rho = 1/2\n', exp(-2*sp^2*(1 - exp(-1))));
fprintf('Patch half-width at S = 1/2 (large aperture): dx = (z/d) l sqrt(-ln(1 - ln2/(2 sig^2))) = %.2f mm at d = 10 mm\n', ...
    z/10e-3*ell*sqrt(-log(1 - log(2)/(2*sp^2)))*1e3);

figure('Position', [50 50 600 380]); cols = 'brk';
hl = zeros(1, numel(d_list) + 1);
for a = 1:numel(d_list)
    hl(a) = plot(dx*1e3, Sex(a, :), [cols(a) '-']); hold on; plot(dx*1e3, Smc(a, :), [cols(a) 'o']);
end
hl(end) = plot(dx*1e3, exp(-sp^2)*ones(size(dx)), 'k:'); xlim([0 16]); xlabel('distance from the guide point (mm)'); ylabel('Strehl S');
title('One law, many points (lines exact, circles MC)'); legend(hl, 'd = 5 mm', 'd = 10 mm', 'd = 20 mm', 'exp(-\sigma^2)');
