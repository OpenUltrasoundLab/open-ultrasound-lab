%% m19_matrix_gratinglobe_uv.m  --  2D grating lobes in (u, v) and the micro-beamformer steering limit
% Purpose : (1) 2D array factor of an Nx x Ny rectangular grid (pitches dx, dy) steered to (th, phi),
%           in direction-cosine space u = sin(th)cos(phi), v = sin(th)sin(phi). Replicas of the main
%           lobe sit at (u0 + m lambda/dx, v0 + n lambda/dy); only those inside the unit circle are
%           real grating lobes. Prints the visible replicas, their level with and without a soft-baffle
%           element factor, and the peak level outside the main lobe for five cases.
%           (2) Largest pitch with no grating lobe for steering anywhere in a cone sin(th) <= s:
%           square grid d <= lambda/(1+s); hexagonal (triangular) grid, nearest-neighbour spacing
%           d <= 2 lambda/(sqrt(3)(1+s)), 13.4 % fewer elements for the same aperture.
%           (3) Micro-beamformer (K x K subarrays at pitch D = K d): closed-form ratio of the worst
%           coarse grating lobe to the main lobe, rho = sin(a)/sin(pi/K - a), a = pi d delta/lambda,
%           and the steering window around the pre-steer for a lobe limit, checked by a direct sum.
% Change  : the case list in (1), cone angles in (2), K list and lobe limits in (3).
% Expect  : (1) 0.5 lambda, (45,0) deg: no visible replica, peak -13.2 dB (main-lobe sidelobe, grid).
%           1.0 lambda, (20,0): replica at u = -0.658, 0 dB isotropic, -6.3 dB with element factor.
%           0.7 lambda, (30,0): replica at u = -0.929, 0 dB / -11.5 dB; (30,45): replicas at radius
%           1.132, outside, peak -13.2 dB. 0.5 x 0.8 lambda, (30,90): replica at v = -0.75, -5.1 dB.
%           (2) cone 30/45/60/90 deg: square 0.667/0.586/0.536/0.500 lambda, hex 0.770/0.676/0.619/0.577.
%           (3) d = lambda/2, 48 x 48: cables 576/256/144/64 for K = 2/3/4/6; lobe at 10 deg off
%           pre-steer -11.1/-8.3/-5.2/+0.7 dB; -20 dB window +-3.6/3.0/2.4/1.7 deg; -10 dB window
%           +-11.2/8.5/6.6/4.5 deg. K = 3 pre-steered to 30 deg: -20 dB beams from 26.6 to 33.5 deg.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.  Runtime about 2 s.
clear; close all; rng(1);

% Dirichlet kernel (array factor of N points at pitch 1 in units of x = d u / lambda), safe at integers
Dn = @(x, N) (abs(sin(pi*x)) > 1e-9) .* sin(N*pi*x) ./ (N*sin(pi*x) + (abs(sin(pi*x)) <= 1e-9)) ...
           + (abs(sin(pi*x)) <= 1e-9) .* cos(N*pi*x) ./ cos(pi*x);
snc = @(x) (abs(x) > 1e-12) .* sin(pi*x) ./ (pi*x + (abs(x) <= 1e-12)) + (abs(x) <= 1e-12);
% soft-baffle rectangular element of width wx, wy (in lambda): sinc(wx u) sinc(wy v) cos(th)
EF = @(u, v, wx, wy) abs(snc(wx*u) .* snc(wy*v)) .* sqrt(max(0, 1 - u.^2 - v.^2));

%% (1) 2D array factor and grating lobes for five cases
% columns: dx/lambda, dy/lambda, th0 (deg), phi0 (deg)
cases = [0.5 0.5 45 0; 1.0 1.0 20 0; 0.7 0.7 30 0; 0.7 0.7 30 45; 0.5 0.8 30 90];
Nx = 32; Ny = 32; du = 0.0025; g = -1:du:1; [U, V] = meshgrid(g, g); vis = U.^2 + V.^2 <= 1;
fprintf('(1) %d x %d grid, uniform weights, soft-baffle element of width 0.9 pitch\n', Nx, Ny);
for c = 1:size(cases, 1)
  dx = cases(c,1); dy = cases(c,2); th = cases(c,3)*pi/180; ph = cases(c,4)*pi/180;
  u0 = sin(th)*cos(ph); v0 = sin(th)*sin(ph);
  AF = abs(Dn(dx*(U - u0), Nx) .* Dn(dy*(V - v0), Ny));
  E = EF(U, V, 0.9*dx, 0.9*dy); Em = EF(u0, v0, 0.9*dx, 0.9*dy);
  main = abs(U - u0) < 1/(Nx*dx) & abs(V - v0) < 1/(Ny*dy);       % first-null box of the main lobe
  pk_iso = 20*log10(max(AF(vis & ~main)));
  pk_ef = 20*log10(max(AF(vis & ~main) .* E(vis & ~main)) / Em);
  fprintf('  dx = %.2f, dy = %.2f lambda, steer (th, phi) = (%g, %g) deg, (u0, v0) = (%.3f, %.3f)\n', ...
          dx, dy, cases(c,3), cases(c,4), u0, v0);
  nvis = 0;
  for m = -3:3
    for n = -3:3
      if m == 0 && n == 0, continue; end
      ur = u0 + m/dx; vr = v0 + n/dy; r = hypot(ur, vr);
      if r <= 1
        nvis = nvis + 1;
        fprintf('     replica (m,n) = (%2d,%2d) at (u, v) = (%6.3f, %6.3f), |.| = %.3f: 0 dB isotropic, %6.1f dB with element factor\n', ...
                m, n, ur, vr, r, 20*log10(EF(ur, vr, 0.9*dx, 0.9*dy)/Em));
      elseif r <= 1.25
        fprintf('     replica (m,n) = (%2d,%2d) at |.| = %.3f: just outside the visible circle\n', m, n, r);
      end
    end
  end
  if nvis == 0, fprintf('     no replica inside the visible circle\n'); end
  fprintf('     peak outside main lobe: %6.1f dB isotropic, %6.1f dB with element factor\n', pk_iso, pk_ef);
end
% check the closed form against a direct phasor sum at random (u, v) for case 4
dx = 0.7; dy = 0.7; u0 = sin(pi/6)*cos(pi/4); v0 = sin(pi/6)*sin(pi/4);
[ix, iy] = meshgrid(0:Nx-1, 0:Ny-1); ix = ix(:); iy = iy(:);
up = 2*rand(50,1) - 1; vp = 2*rand(50,1) - 1;
Pd = abs(exp(1i*2*pi*(dx*(up - u0)*ix.' + dy*(vp - v0)*iy.')) * ones(Nx*Ny,1)) / (Nx*Ny);
Pc = abs(Dn(dx*(up - u0), Nx) .* Dn(dy*(vp - v0), Ny));
fprintf('  closed form vs direct %d-element sum: max error %.1e\n', Nx*Ny, max(abs(Pd - Pc)));

%% (2) Largest pitch with no grating lobe over a steering cone
fprintf('\n(2) Max pitch for no grating lobe, steering anywhere with th <= th_max\n');
fprintf('   th_max   square d/lambda   hex d/lambda (nearest neighbour)\n');
for tm = [30 45 60 90]
  s = sin(tm*pi/180);
  fprintf('   %4d     %.3f             %.3f\n', tm, 1/(1+s), 2/(sqrt(3)*(1+s)));
end
% widest cone for a given pitch: sin(th_max) = |g_min| - 1 (|g_min| = lambda/d square, 2 lambda/(sqrt(3) d) hex)
cone = @(gmin) asin(min(1, max(0, gmin - 1)))*180/pi;
fprintf('   d = 0.70 lambda: no-GL cone %.1f deg (square), %.1f deg (hex, same nearest-neighbour spacing)\n', cone(1/0.7), cone(2/(sqrt(3)*0.7)));
dh = 0.6*sqrt(2/sqrt(3));                                 % hex spacing with the element density of a 0.6 lambda square grid
fprintf('   d = 0.60 lambda square: cone %.1f deg; hex at equal density (d = %.3f lambda): cone %.1f deg\n', cone(1/0.6), dh, cone(2/(sqrt(3)*dh)));
% reciprocal lattice of the triangular grid, a1 = d(1,0), a2 = d(1/2, sqrt(3)/2), d = 1 lambda
a1 = [1 0]; a2 = [0.5 sqrt(3)/2]; B = inv([a1; a2]).';          % rows b1, b2 with a_i . b_j = delta_ij
best = inf; for m = -2:2, for n = -2:2, if m || n, best = min(best, norm(m*B(1,:) + n*B(2,:))); end, end, end
hx = [a1; a2]; pts = []; for i = -6:6, for j = -6:6, pts = [pts; i*a1 + j*a2]; end, end %#ok<AGROW>
AFhex = abs(sum(exp(1i*2*pi*(pts*B(1,:).')))) / size(pts,1);
fprintf('   triangular grid: nearest replica at %.4f lambda/d (2/sqrt(3) = %.4f); |AF| there = %.3f\n', best, 2/sqrt(3), AFhex);
fprintf('   element density hex/square at equal no-GL cone: %.3f (%.1f %% fewer)\n', sqrt(3)/2, 100*(1 - sqrt(3)/2));

%% (3) Micro-beamformer: worst coarse grating lobe vs pre-steer error, and the steering window
dl = 0.5; N = 48;                                       % element pitch d/lambda, elements per axis
rho = @(dlt, K) sin(pi*dl*dlt) ./ sin(pi/K - pi*dl*dlt); % lobe/main, delta in u units (0 <= delta < 1/(K dl))
win = @(r, K) atan(r*sin(pi/K) ./ (1 + r*cos(pi/K))) / (pi*dl);   % delta_max in u for a lobe ratio r
fprintf('\n(3) Micro-beamforming, d = lambda/2, %d x %d matrix, pre-steer broadside\n', N, N);
fprintf('   K   cables   D/lambda  coarse GL spacing   lobe at 10 deg   window -10 dB   window -20 dB\n');
u = linspace(-1, 1, 8001); Kl = [2 3 4 6];
for K = Kl
  dlt = sin(10*pi/180);
  w10 = asin(win(10^(-10/20), K))*180/pi; w20 = asin(win(10^(-20/20), K))*180/pi;
  fprintf('   %d   %4d     %.1f       %.3f             %6.1f dB        +-%.1f deg       +-%.1f deg\n', ...
          K, (N/K)^2, K*dl, 1/(K*dl), 20*log10(rho(dlt, K)), w10, w20);
end
% direct-sum check (1D cut): K = 3, pre-steer 0, beam 10 deg
K = 3; n = 0:N-1; xi = mod(n, K) - (K-1)/2; Xk = (floor(n/K) + 0.5)*K - N/2;
us = sin(10*pi/180); w = exp(-1i*2*pi*dl*Xk*us);              % fine delays: pre-steer 0 (none)
P = abs(exp(1i*2*pi*dl*u(:)*(Xk + xi)) * w(:)).'; P = P/max(P);
gl = max(P(abs(u - (us - 1/(K*dl))) < 0.05));
fprintf('   direct sum, K = 3, beam 10 deg: lobe at u = %.3f is %.2f dB (closed form %.2f dB)\n', ...
        us - 1/(K*dl), 20*log10(gl), 20*log10(rho(us, K)));
% pre-steer at 30 deg: the window is symmetric in u, not in angle
K = 3; u0 = 0.5; dm = win(0.1, K);
fprintf('   K = 3, pre-steer 30 deg, -20 dB: beams from %.1f to %.1f deg (u %.3f to %.3f)\n', ...
        asin(u0 - dm)*180/pi, asin(u0 + dm)*180/pi, u0 - dm, u0 + dm);
% 2D: steering off in both axes; the diagonal coarse replica scales as rho_x * rho_y
dlt = sin(5*pi/180);
fprintf('   2D, K = 3, delta = 5 deg in u and in v: axis lobes %.1f dB, diagonal lobe %.1f dB\n', ...
        20*log10(rho(dlt, 3)), 40*log10(rho(dlt, 3)));
% number of pre-steer zones to cover a +-45 deg sector at -20 dB (zones tile u uniformly)
for K = [2 3 4]
  fprintf('   K = %d: %d pre-steer zones per axis cover u in [-0.707, 0.707] at -20 dB\n', K, ceil(2*0.7071/(2*win(0.1, K))));
end

%% Figures
figure('Position', [50 50 1000 420]);
subplot(1, 2, 1);
gg = -2:0.01:2; [UU, VV] = meshgrid(gg, gg); dx = 0.7; dy = 0.7; u0 = 0.5; v0 = 0;
A = abs(Dn(dx*(UU - u0), 16) .* Dn(dy*(VV - v0), 16));
imagesc(gg, gg, max(-40, 20*log10(A + 1e-12))); axis image; set(gca, 'YDir', 'normal'); colorbar; hold on;
t = linspace(0, 2*pi, 200); plot(cos(t), sin(t), 'w-', 'LineWidth', 1.5);
[mm, nn] = meshgrid(-3:3, -3:3); plot(u0 + mm(:)/dx, v0 + nn(:)/dy, 'w+');
xlim([-2 2]); ylim([-2 2]); xlabel('u = sin\theta cos\phi'); ylabel('v = sin\theta sin\phi');
title('16 x 16 at 0.7\lambda, steer (30^o, 0^o): AF (dB), visible circle, replicas');
subplot(1, 2, 2);
cols = 'brkm'; hold on; hl = zeros(1, numel(Kl));
for i = 1:numel(Kl)
  K = Kl(i); th = linspace(0, asin(min(0.999, 0.98/(K*dl)))*180/pi, 200);
  hl(i) = plot(th, 20*log10(rho(sin(th*pi/180), K)), cols(i), 'LineWidth', 1.5);
end
plot([0 30], [-20 -20], 'k:'); plot([0 30], [-10 -10], 'k:');
xlim([0 30]); ylim([-40 5]); grid on; xlabel('beam minus pre-steer (deg), pre-steer broadside');
ylabel('worst coarse grating lobe / main (dB)'); legend(hl, 'K = 2', 'K = 3', 'K = 4', 'K = 6', 'Location', 'southeast');
title('Micro-beamformer steering limit, d = \lambda/2');
