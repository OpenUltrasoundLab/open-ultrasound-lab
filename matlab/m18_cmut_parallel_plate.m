%% m18_cmut_parallel_plate.m
% CMUT electromechanics: lumped parallel-plate model versus a distributed
% clamped-plate model (axisymmetric finite differences), conventional and
% collapse branches, snap-back, spring softening and coupling versus bias.
%
% What it shows
%   1) Lumped piston: (V/Vpi)^2 = 27/4 xi (1-xi)^2, kT^2 = 2xi/(1-xi), pull-in at xi = 1/3.
%   2) Distributed clamped plate D*del^4 w = eps0 V^2 / (2 (g_eff - w)^2):
%      conventional branch traced by fixing the centre deflection (the voltage
%      that holds it is the unknown), collapsed branch traced by fixing the
%      contact radius (zero contact reaction at its edge). Both are Newton solves.
%   3) Pull-in voltage, centre and average deflection at pull-in, snap-back voltage.
%   4) Coupling k^2 = 1 - lambda_V / lambda_Q from the two lowest eigenvalues
%      (constant voltage vs constant charge) on the conventional branch.
% What to change: a, h, g0, ti, er (geometry), N (grid), E/nu/rho (membrane).
% Expected output (defaults: Si membrane a = 20 um, h = 1 um, g0 = 0.2 um,
%   0.1 um nitride): f0 = 10.24 MHz, Vpi = 53.1 V distributed vs 76.0 V lumped
%   (ratio 0.699), centre deflection at pull-in 0.46 g_eff (0.459 printed below), average 0.15 g_eff,
%   snap-back about 26 V, k^2 = 0.26 at 0.9 Vpi (lumped model: 0.43).
% Runtime: about 5 s in Octave.
clear; close all; rng(1);

%% Parameters (SI units)
E = 169e9; nu = 0.22; rho = 2330;     % silicon membrane
a = 20e-6; h = 1.0e-6;                % membrane radius, thickness (m)
g0 = 0.20e-6; ti = 0.10e-6; er = 7.5; % vacuum gap, insulation thickness, insulation permittivity
e0 = 8.854e-12;
N = 60;                               % radial grid intervals
D = E*h^3/(12*(1-nu^2));              % flexural rigidity (N m)
ge = g0 + ti/er;                      % effective electrical gap (m)
f0 = 10.2158/(2*pi*a^2)*sqrt(D/(rho*h));
fprintf('Membrane: a = %.0f um, h = %.2f um, D = %.3g N m, f0 (vacuum) = %.2f MHz\n', a*1e6, h*1e6, D, f0/1e6);

%% Lumped piston model (spring from the average deflection under uniform pressure)
kp = 192*pi*D/a^2; Ap = pi*a^2;
Vpi_l = sqrt(8*kp*ge^3/(27*e0*Ap));
xi = linspace(0, 1/3, 400);
Vl = sqrt(27/4*xi.*(1-xi).^2);        % V/Vpi
k2l = 2*xi./(1-xi); fl = sqrt((1-3*xi)./(1-xi));
fprintf('Lumped: Vpi = %.1f V; at 0.9 Vpi: xi = %.3f, kT^2 = %.3f, f/f0 = %.3f\n', Vpi_l, ...
  interp1(Vl, xi, 0.9), interp1(Vl, k2l, 0.9), interp1(Vl, fl, 0.9));

%% Distributed model: axisymmetric FD operators
dr = a/N; r = (0:N)'*dr;
L = zeros(N+1, N);                    % Laplacian at nodes 0..N acting on unknowns w(0..N-1); w(N) = 0, w'(a) = 0
L(1,1) = -4/dr^2; L(1,2) = 4/dr^2;
for i = 2:N
  rp = r(i) + dr/2; rm = r(i) - dr/2;
  L(i,i) = -(rp+rm)/(r(i)*dr^2); L(i,i-1) = rm/(r(i)*dr^2);
  if i < N, L(i,i+1) = rp/(r(i)*dr^2); end
end
L(N+1, N) = 2/dr^2;                   % clamped edge via ghost node
Wt = 2*pi*r*dr; Wt(1) = pi*(dr/2)^2; Wt(end) = pi*r(end)*dr;
K = D*(L'*diag(Wt)*L); Wf = Wt(1:N); M = rho*h*diag(Wf);
lam0 = min(real(eig(K, M)));
fprintf('FD check: f0 = %.3f MHz (exact 10.2158 constant gives %.3f MHz)\n', sqrt(lam0)/2/pi/1e6, f0/1e6);

modal = @(KV, KQ, F) [min(real(eig(KV(F,F), M(F,F)))), min(real(eig(KQ(F,F), M(F,F))))];

%% Conventional branch: control the centre deflection wc, unknowns w and lambda = V^2
wcs = linspace(0.01, 0.99, 99)*g0;
w = zeros(N,1); lam = 0;
Vc = zeros(size(wcs)); k2c = Vc; fc = Vc; wav = Vc;
for s = 1:numel(wcs)
  for it = 1:50
    fe = Wf*e0./(2*(ge-w).^2); fp = Wf*e0./(ge-w).^3;
    R = [K*w - lam*fe; w(1) - wcs(s)];
    J = [K - lam*diag(fp), -fe; [1, zeros(1,N-1)], 0];
    dx = -J\R; w = w + dx(1:N); lam = lam + dx(end);
    if max(abs(dx(1:N))) < 1e-15, break; end
  end
  Vc(s) = sqrt(lam); wav(s) = sum(Wf.*w)/sum(Wt);
  fp = lam*Wf*e0./(ge-w).^3; KV = K - diag(fp);
  gC = Wf*e0./(ge-w).^2; C = sum(Wf*e0./(ge-w));
  KQ = KV + lam*(gC*gC')/C;
  ev = modal(KV, KQ, 1:N);
  k2c(s) = 1 - ev(1)/ev(2); fc(s) = sqrt(max(ev(1),0))/2/pi;
  if s == 1, wsave = zeros(N,0); end
  if any(abs(wcs(s)/g0 - [0.2 0.4 0.49]) < 1e-9), wsave(:,end+1) = w; end %#ok<AGROW>
end
[Vpi, ip] = max(Vc);
fprintf('Distributed: Vpi = %.1f V (lumped %.1f V, ratio %.3f)\n', Vpi, Vpi_l, Vpi/Vpi_l);
fprintf('  at pull-in: centre deflection %.3f g_eff, average deflection %.3f g_eff\n', wcs(ip)/ge, wav(ip)/ge);
fprintf('  closed form fit: Vpi^2 = %.2f D g_eff^3 / (eps0 a^4)\n', Vpi^2*e0*a^4/(D*ge^3));
st = 1:ip;                             % stable part of the conventional branch
for fr = [0.5 0.8 0.9 0.95]
  fprintf('  V/Vpi = %.2f: k^2 = %.3f, f/f0 = %.3f, wc/g0 = %.3f\n', fr, interp1(Vc(st)/Vpi, k2c(st), fr), ...
    interp1(Vc(st)/Vpi, fc(st)/f0, fr), interp1(Vc(st)/Vpi, wcs(st)/g0, fr));
end

%% Collapsed branch: control the contact radius (nodes 1..nc at w = g0, zero reaction at node nc)
w(1) = g0; Vcol = []; bcol = []; wavc = [];
for nc = 1:N-4
  F = (nc+1:N)'; e = nc;
  for it = 1:60
    w(1:nc) = g0;
    fe = Wf*e0./(2*(ge-w).^2); fp = Wf*e0./(ge-w).^3;
    Kf = K - lam*diag(fp); Rv = K*w - lam*fe;
    R = [Rv(F); Rv(e)];
    J = [Kf(F,F), -fe(F); Kf(e,F), -fe(e)];
    dx = -J\R; sstep = 1;
    while any(w(F) + sstep*dx(1:end-1) >= ge - 1e-12) && sstep > 1e-4, sstep = sstep/2; end
    w(F) = w(F) + sstep*dx(1:end-1); lam = lam + sstep*dx(end);
    if max(abs(sstep*dx(1:end-1))) < 1e-15, break; end
  end
  Vcol(end+1) = sqrt(lam); bcol(end+1) = r(nc); wavc(end+1) = sum(Wf.*w)/sum(Wt); %#ok<AGROW>
end
[Vsb, isb] = min(Vcol);
fprintf('Collapsed branch: snap-back at %.1f V (%.3f Vpi); contact radius b/a = 0.5 at %.2f Vpi\n', ...
  Vsb, Vsb/Vpi, interp1(bcol(isb:end)/a, Vcol(isb:end), 0.5)/Vpi);
fprintf('Hysteresis window: %.1f V to %.1f V\n', Vsb, Vpi);

%% Figures
figure('Position', [80 80 1100 380]);
subplot(1,3,1);
plot(Vc(st), wav(st)/ge, 'b-', 'LineWidth', 2); hold on;
plot(Vc(ip:end), wav(ip:end)/ge, 'b--');
plot(Vcol, wavc/ge, 'r-', 'LineWidth', 2);
plot(Vl*Vpi_l, xi, 'k:', 'LineWidth', 1.5);
xlim([0 2*Vpi]); ylim([0 1]);
xlabel('bias V (V)'); ylabel('average deflection / g_{eff}');
title('Both branches and hysteresis'); legend('conventional, stable', 'conventional, unstable', 'collapsed', 'lumped piston', 'Location', 'northwest');
subplot(1,3,2);
plot(Vl, k2l, 'k:', 'LineWidth', 1.5); hold on; plot(Vc(st)/Vpi, k2c(st), 'b-', 'LineWidth', 2);
xlabel('V / V_{PI}'); ylabel('k^2'); ylim([0 1]); title('Coupling vs bias'); legend('lumped piston', 'distributed', 'Location', 'northwest');
subplot(1,3,3);
plot(r(1:N)/a*1e0, wsave/g0, 'LineWidth', 1.5); xlabel('r / a'); ylabel('w / g_0'); title('Profiles at w_c = 0.2, 0.4, 0.49 g_0');
