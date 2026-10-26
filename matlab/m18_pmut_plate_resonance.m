%% m18_pmut_plate_resonance.m
% PMUT: clamped circular plate resonance, bilayer (piezo film on passive layer),
% water loading (Lamb added mass) and the optimum inner-electrode radius.
%
% What it shows
%   1) The clamped-plate eigenvalue: J0(L) I1(L) + I0(L) J1(L) = 0 gives L = 3.1962,
%      L^2 = 10.2158; an FD eigen-solve of the same plate converges to it.
%   2) f = (L^2 / (2 pi a^2)) sqrt(D / (rho h)) for a bilayer with the neutral-axis D.
%   3) Water on one side: f_w = f_vac / sqrt(1 + 0.67 rho_w a / (rho h))  (Lamb 1920).
%   4) Actuation from a uniform piezoelectric moment over r < beta*a is proportional to
%      beta^2 (1 - beta^2) (one-term Ritz): optimum beta = 0.707, zero for a full electrode.
% What to change: film material, film and passive thickness, radius list.
% Expected output: lambda^2 = 10.2158 (FD 10.2113, 10.2147, 10.2155 for N = 40, 80, 160);
%   1 um film on 2 um Si, a = 50 um: f_vac = 5.58 / 4.81 / 2.98 MHz (AlN / ScAlN / PZT),
%   in water 2.44 / 2.11 / 1.54 MHz (Lamb factor 4.2, 4.2, 2.7); optimum beta = 0.707.
% Runtime: about 2 s.
clear; close all; rng(1);

%% 1) Clamped plate eigenvalue from the Bessel characteristic equation
g = @(L) besselj(0,L).*besseli(1,L) + besseli(0,L).*besselj(1,L);
lo = 3.0; hi = 3.4;                     % bracket around the first root
for it = 1:60, mid = (lo+hi)/2; if sign(g(mid)) == sign(g(lo)), lo = mid; else hi = mid; end, end
L1 = (lo+hi)/2;
fprintf('Clamped circular plate: lambda = %.4f, lambda^2 = %.4f\n', L1, L1^2);

%% 2) FD check of the constant (unit plate, D = rho h = 1, a = 1)
for N = [40 80 160]
  dr = 1/N; r = (0:N)'*dr; Lp = zeros(N+1,N);
  Lp(1,1) = -4/dr^2; Lp(1,2) = 4/dr^2;
  for i = 2:N
    rp = r(i)+dr/2; rm = r(i)-dr/2;
    Lp(i,i) = -(rp+rm)/(r(i)*dr^2); Lp(i,i-1) = rm/(r(i)*dr^2);
    if i < N, Lp(i,i+1) = rp/(r(i)*dr^2); end
  end
  Lp(N+1,N) = 2/dr^2;
  Wt = 2*pi*r*dr; Wt(1) = pi*(dr/2)^2; Wt(end) = pi*dr;
  K = Lp'*diag(Wt)*Lp; M = diag(Wt(1:N));
  fprintf('  FD N = %3d: lambda^2 = %.4f\n', N, sqrt(min(real(eig(K,M)))));
end

%% 3) Bilayer plates: film on silicon
% Film data: e31,f and eps_r from Wang et al. (Hilton Head 2016) / Roy et al. (2023);
% E as listed by Wang et al.; densities are typical bulk values.
films = {'AlN', 390e9, 0.24, 3260, -1.05, 10.5;  'ScAlN', 200e9, 0.24, 3300, -1.6, 12; ...
         'PZT film', 65e9, 0.30, 7600, -14.0, 1200};
Es = 169e9; nus = 0.22; rhos = 2330;    % silicon passive layer
rhow = 1000;
hs = 2e-6; hp = 1e-6; a = 50e-6;
fprintf('\nBilayer: %.1f um film on %.1f um Si, a = %.0f um\n', hp*1e6, hs*1e6, a*1e6);
fprintf('%-9s %8s %8s %9s %9s %10s %12s\n', 'film', 'zn (um)', 'zp (um)', 'fvac MHz', 'fwat MHz', 'beta_Lamb', 'wc/V (nm/V)');
res = zeros(size(films,1), 3);
for k = 1:size(films,1)
  Ep = films{k,2}; nup = films{k,3}; rhop = films{k,4}; e31 = films{k,5};
  E1 = Es/(1-nus^2); E2 = Ep/(1-nup^2);          % plate moduli
  z1 = hs/2; z2 = hs + hp/2;                      % layer mid-planes from the bottom
  zn = (E1*hs*z1 + E2*hp*z2)/(E1*hs + E2*hp);     % neutral axis
  Db = E1*(hs^3/12 + hs*(z1-zn)^2) + E2*(hp^3/12 + hp*(z2-zn)^2);
  mh = rhos*hs + rhop*hp;                         % mass per area
  fv = L1^2/(2*pi*a^2)*sqrt(Db/mh);
  bL = 0.67*rhow*a/mh; fw = fv/sqrt(1 + bL);
  zp = z2 - zn;                                   % moment arm of the film
  wcV = 0.375*abs(e31)*zp*a^2*0.25/Db;            % centre deflection per volt at beta^2 = 1/2 (one-term Ritz)
  fprintf('%-9s %8.2f %8.2f %9.2f %9.2f %10.2f %12.2f\n', films{k,1}, zn*1e6, zp*1e6, fv/1e6, fw/1e6, bL, wcV*1e9);
  res(k,:) = [fv fw wcV];
end

%% 4) Electrode coverage factor
beta = linspace(0,1,2001); q = beta.^2.*(1-beta.^2);
[~, ib] = max(q); fprintf('\nOptimum inner-electrode radius: beta = %.3f a (factor %.3f); full electrode gives %.3f\n', beta(ib), q(ib), q(end));

%% 5) Frequency versus radius for the three films, in water
aa = linspace(20e-6, 200e-6, 100); F = zeros(size(films,1), numel(aa));
for k = 1:size(films,1)
  Ep = films{k,2}; nup = films{k,3}; rhop = films{k,4};
  E1 = Es/(1-nus^2); E2 = Ep/(1-nup^2); z1 = hs/2; z2 = hs+hp/2;
  zn = (E1*hs*z1 + E2*hp*z2)/(E1*hs + E2*hp);
  Db = E1*(hs^3/12 + hs*(z1-zn)^2) + E2*(hp^3/12 + hp*(z2-zn)^2); mh = rhos*hs + rhop*hp;
  F(k,:) = L1^2./(2*pi*aa.^2).*sqrt(Db/mh)./sqrt(1 + 0.67*rhow*aa/mh);
end
fprintf('Water-loaded f at a = 100 um: AlN %.2f, ScAlN %.2f, PZT %.2f MHz\n', interp1(aa,F(1,:),100e-6)/1e6, ...
  interp1(aa,F(2,:),100e-6)/1e6, interp1(aa,F(3,:),100e-6)/1e6);

figure('Position', [80 80 1000 380]);
subplot(1,2,1); loglog(aa*1e6, F/1e6, 'LineWidth', 1.8); grid on;
xlabel('radius a (\mum)'); ylabel('f in water (MHz)'); legend(films(:,1)); title('PMUT resonance, 2 \mum Si + 1 \mum film');
subplot(1,2,2); plot(beta, q/max(q), 'LineWidth', 2); grid on; xlabel('electrode radius / a'); ylabel('relative actuation');
title('Inner electrode: optimum at 0.707 a');
