%% m18_composite_effective_properties.m
% 1-3 piezocomposites: Smith-Auld effective thickness-mode properties for
% PZT-5H pillars and PMN-28PT single-crystal pillars in the same soft epoxy.
%
% What it shows: kt, Z, v and clamped permittivity versus ceramic volume fraction,
%   the plateau of kt near nu = 0.4-0.7, and how much a crystal pillar adds.
% What to change: polymer (vL, vS, rho), pillar constants, f0 for the thickness table.
% Constants: PZT-5H as in Module 16 (c11E 12.6, c12E 7.95, c13E 8.41, c33E 11.7 x1e10 N/m^2,
%   e31 -6.5, e33 23.3 C/m^2, eps33S 1470, rho 7500); PMN-28PT from Kim, Joh & Roh,
%   Sensors and Materials 25 (2013) 539 (c11E 14.7, c12E 13.1, c13E 11.4, c33E 11.4,
%   e31 -6.88, e33 21.4, eps33S 885, rho 8084). Soft epoxy vL 2500, vS 1100 m/s, rho 1150.
% Expected output: nu = 1: kt 0.513 / 0.582, Z 34.5 / 37.3 MRayl;
%   nu = 0.5: kt 0.687 / 0.833, Z 16.4 / 17.0 MRayl.
% Runtime: < 1 s.
clear; close all; rng(1);
e0 = 8.854e-12;
pz = [12.6e10 7.95e10 8.41e10 11.7e10 -6.5 23.3 1470 7500];   % c11 c12 c13 c33 e31 e33 epsS rho
pm = [14.7e10 13.1e10 11.4e10 11.4e10 -6.88 21.4 885 8084];
rp = 1150; vL = 2500; vS = 1100; er_p = 4;
p11 = rp*vL^2; p44 = rp*vS^2; p12 = p11 - 2*p44;

nus = linspace(0.05, 1, 96);
out = zeros(2, numel(nus), 4);
mats = {pz, pm};
for m = 1:2
  c = mats{m};
  for i = 1:numel(nus)
    nu = nus(i);
    Dn  = nu*(p11 + p12) + (1-nu)*(c(1) + c(2));
    c33E = nu*(c(4) - 2*(1-nu)*(c(3) - p12)^2/Dn) + (1-nu)*p11;
    e33 = nu*(c(6) - 2*(1-nu)*c(5)*(c(3) - p12)/Dn);
    eS  = nu*(c(7)*e0 + 2*(1-nu)*c(5)^2/Dn) + (1-nu)*er_p*e0;
    c33D = c33E + e33^2/eS; rho = nu*c(8) + (1-nu)*rp;
    v = sqrt(c33D/rho);
    out(m,i,:) = [sqrt(e33^2/(c33D*eS)), rho*v/1e6, v, eS/e0];
  end
end
fprintf('%5s | %18s | %18s\n', 'nu', 'PZT-5H: kt  Z  epsS', 'PMN-28PT: kt  Z  epsS');
for nu = [1 0.7 0.6 0.5 0.4 0.3 0.2]
  [~, i] = min(abs(nus - nu));
  fprintf('%5.2f | %5.3f %5.1f %6.0f | %5.3f %5.1f %6.0f\n', nus(i), out(1,i,1), out(1,i,2), out(1,i,4), out(2,i,1), out(2,i,2), out(2,i,4));
end
f0 = 5e6; [~, i] = min(abs(nus - 0.5));
fprintf('At nu = 0.5, f0 = 5 MHz: half-wave thickness %.0f um (PZT) and %.0f um (PMN-PT); pillar width <= 0.5 t: %.0f um\n', ...
  out(1,i,3)/(2*f0)*1e6, out(2,i,3)/(2*f0)*1e6, 0.5*out(2,i,3)/(2*f0)*1e6);

figure('Position', [80 80 1000 380]);
subplot(1,2,1); plot(nus, squeeze(out(1,:,1)), 'b-', nus, squeeze(out(2,:,1)), 'r-', 'LineWidth', 2); grid on;
xlabel('ceramic volume fraction \nu'); ylabel('k_t'); legend('PZT-5H / epoxy', 'PMN-28PT / epoxy', 'Location', 'south'); ylim([0.3 0.9]);
title('Effective thickness coupling (Smith-Auld)');
subplot(1,2,2); plot(nus, squeeze(out(1,:,2)), 'b-', nus, squeeze(out(2,:,2)), 'r-', 'LineWidth', 2); grid on;
xlabel('ceramic volume fraction \nu'); ylabel('Z (MRayl)'); title('Acoustic impedance');
