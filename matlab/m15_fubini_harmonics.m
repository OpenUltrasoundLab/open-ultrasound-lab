%% m15_fubini_harmonics.m
% Purpose : plane-wave nonlinear propagation. Shock distance xbar = rho c^3/(beta omega p0),
%           Fubini harmonics B_n = 2 J_n(n sigma)/(n sigma), and a coupled-harmonic solution of
%           the Burgers equation with power-law absorption (Haran & Cook style):
%             dB_n/dsigma = (n/4) sum_{k=1}^{n-1} B_k B_{n-k} - (n/2) sum_{k>=1} B_k B_{k+n} - alpha_n xbar B_n
%           where p/p0 = sum B_n sin(n omega tau), sigma = x/xbar.
% Change  : medium table, p0, f0, number of harmonics N.
% Expect  : water (B/A 5.0) 3 MHz 1 MPa: xbar = 49.2 mm; liver (B/A 7.0) 49.1 mm; fat (B/A 10) 25.6 mm
%           Fubini at sigma = 0.5: B2 = 0.230 (B2/B1 = -12.5 dB); ODE lossless matches to 1e-4
%           liver with absorption (Gamma = 1.06): B2 peaks near 24 mm (23.6) at -20.9 dB re p0
% Runtime : about 3 s.
clear; close all; rng(1);

%% Media: rho (kg/m^3), c (m/s), alpha0 (dB/cm/MHz^y), y, B/A
names = {'water 20 C', 'blood', 'liver', 'muscle', 'fat'};
M = [ 998 1482 0.0022 2.0  5.0;
     1060 1584 0.15   1.2  6.1;
     1060 1578 0.5    1.1  7.0;
     1050 1580 1.0    1.1  7.4;
      950 1450 0.6    1.0 10.0];
f0 = 3e6; p0 = 1e6; dBperNp = 20*log10(exp(1));
fprintf('%-11s  beta   xbar(mm)  Mach eps   alpha(f0) Np/m   Gamma = 1/(alpha xbar)\n', 'medium');
for i = 1:numel(names)
  rho = M(i,1); c = M(i,2); beta = 1 + M(i,5)/2;
  xbar = rho*c^3/(beta*2*pi*f0*p0);
  alpha = M(i,3)*(f0/1e6)^M(i,4)/dBperNp*100;
  fprintf('%-11s  %4.2f   %6.1f    %.2e      %7.2f          %6.2f\n', names{i}, beta, xbar*1e3, p0/(rho*c)/c, alpha, 1/(alpha*xbar));
end

%% Fubini versus the coupled-harmonic ODE (lossless), RK4
N = 40;                                  % harmonics kept
nn = (1:N)';
fub = @(n, s) 2*besselj(n, n*s)./(n*s);
s1 = @(cc) [0; cc(1:N-1)];               % sum_{k=1}^{n-1} B_k B_{n-k}  from conv(B,B)
s2 = @(r) [r(N-1:-1:1); 0];              % sum_k B_k B_{k+n}           from conv(B,flipud(B))
dBds = @(B, att) (nn/4).*s1(conv(B, B)) - (nn/2).*s2(conv(B, flipud(B))) - att.*B;
B = [1; zeros(N-1,1)]; h = 0.0025; att0 = zeros(N,1); sig = 0;
for k = 1:round(0.9/h)
  k1 = dBds(B, att0); k2 = dBds(B + h/2*k1, att0); k3 = dBds(B + h/2*k2, att0); k4 = dBds(B + h*k3, att0);
  B = B + h/6*(k1 + 2*k2 + 2*k3 + k4); sig = sig + h;
  if abs(sig - 0.5) < h/2, B05 = B; end
end
fprintf('\nLossless check (N = %d):   B1       B2       B3\n', N);
fprintf(' Fubini  sigma = 0.5   %.4f   %.4f   %.4f\n', fub(1,0.5), fub(2,0.5), fub(3,0.5));
fprintf(' ODE     sigma = 0.5   %.4f   %.4f   %.4f\n', B05(1), B05(2), B05(3));
fprintf(' Fubini  sigma = 0.9   %.4f   %.4f   %.4f\n', fub(1,0.9), fub(2,0.9), fub(3,0.9));
fprintf(' ODE     sigma = 0.9   %.4f   %.4f   %.4f\n', B(1), B(2), B(3));
fprintf(' small sigma: B2 ~ sigma/2 ; at sigma = 0.1: Fubini %.4f\n', fub(2, 0.1));

%% Liver with power-law absorption: split step (exact absorption, RK4 nonlinearity)
N = 60; nn = (1:N)';
s1 = @(cc) [0; cc(1:N-1)]; s2 = @(r) [r(N-1:-1:1); 0];
dBnl = @(B) (nn/4).*s1(conv(B, B)) - (nn/2).*s2(conv(B, flipud(B)));
rho = 1060; c = 1578; a0 = 0.5; y = 1.1; beta = 1 + 7.0/2;
xbar = rho*c^3/(beta*2*pi*f0*p0);
alphan = a0*(nn*f0/1e6).^y/dBperNp*100;          % Np/m for each harmonic
xmax = 0.12; h = 0.005; ns = round(xmax/xbar/h); hx = xmax/ns; hs = hx/xbar;
E = exp(-alphan*hx/2);                            % half-step absorption factor
B = [1; zeros(N-1,1)]; X = (0:ns)*hx; H = zeros(4, ns+1); H(:,1) = B(1:4);
for k = 1:ns
  B = E.*B;
  k1 = dBnl(B); k2 = dBnl(B + hs/2*k1); k3 = dBnl(B + hs/2*k2); k4 = dBnl(B + hs*k3);
  B = B + hs/6*(k1 + 2*k2 + 2*k3 + k4);
  B = E.*B; H(:,k+1) = B(1:4);
end
[b2max, im] = max(abs(H(2,:)));
fprintf('\nLiver, 3 MHz, 1 MPa: xbar = %.1f mm, Gamma = %.2f\n', xbar*1e3, 1/(a0*3^y/dBperNp*100*xbar));
fprintf(' x (mm)   B1 (dB re p0)   B2 (dB re p0)   B2/B1 (dB)\n');
for xm = [10 25 50 75 100]
  [~, j] = min(abs(X - xm/1e3));
  fprintf(' %5.0f      %6.1f          %6.1f         %6.1f\n', xm, 20*log10(abs(H(1,j))), 20*log10(abs(H(2,j))), 20*log10(abs(H(2,j)/H(1,j))));
end
fprintf(' B2 maximum %.1f dB re p0 at x = %.1f mm\n', 20*log10(b2max), X(im)*1e3);

%% Figures
figure('Position', [60 60 900 600]);
subplot(2,2,1); ss = linspace(1e-3, 1, 200); hold on;
for n = 1:4, plot(ss, fub(n, ss), 'LineWidth', 1.3); end
xlabel('\sigma = x / xbar'); ylabel('B_n'); legend('n=1','n=2','n=3','n=4'); grid on; title('Fubini (lossless)');
subplot(2,2,2); plot(X*1e3, 20*log10(abs(H')), 'LineWidth', 1.3); ylim([-60 0]);
xlabel('x (mm)'); ylabel('dB re p_0'); legend('f_0','2f_0','3f_0','4f_0'); grid on; title('Liver, 3 MHz, 1 MPa, with absorption');
subplot(2,1,2); th = linspace(0, 4*pi, 800); hold on;
for sg = [0 0.5 0.95]
  % implicit Poisson solution P = sin(theta + sigma P), by fixed-point iteration (converges for sigma < 1)
  P = sin(th); for it = 1:200, P = sin(th + sg*P); end
  plot(th/(2*pi), P, 'LineWidth', 1.3);
end
xlabel('\omega\tau / 2\pi (periods)'); ylabel('p / p_0'); legend('\sigma = 0', '0.5', '0.95 (near shock)'); grid on;
title('Waveform steepening (lossless, implicit solution)');
