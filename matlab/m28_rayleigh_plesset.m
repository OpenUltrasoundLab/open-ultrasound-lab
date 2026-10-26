%% m28_rayleigh_plesset.m  --  Rayleigh-Plesset (Marmottant form) with a hand-written RK4
% Purpose : integrate the radial dynamics of a free gas bubble and of a lipid-coated
%           bubble (Marmottant buckling/rupture shell), compute the emitted pressure
%           p_s = rho (R^2 R'' + 2 R R'^2) / r and its spectrum, and read off the
%           harmonic, subharmonic and compression-only signatures.
% Model   : rho (R R'' + 3/2 R'^2) = Pg0 (R0/R)^(3g) (1 - 3g R'/c) - P0 - 2 sig(R)/R
%                                   - 4 mu R'/R - 4 ks R'/R^2 - Pac(t)
%           sig(R) = 0 (R < Rb), chi (R^2/Rb^2 - 1) (Rb <= R <= Rr), sigma_w (R > Rr).
% Change  : R0, f, P (cases table below), chi, ks, s0 (surface tension at rest).
% Expect  : free 2 um bubble (f0 = 2.04 MHz) at 2 MHz, 200 kPa: Rmax/R0 = 2.30,
%           2nd harmonic -10.6 dB; driven at 2 f0 (4.08 MHz) the subharmonic is
%           -30 dB at 100 kPa and -4 dB at 200 kPa (above threshold); a buckled
%           shell (s0 = 0) compresses 2.5 times further than it expands.
% Provenance: simulation (this script). Shell values chi = 0.5 N/m, s0 = 0.02 N/m, ks = 5e-9 kg/s,
%           g = 1.07 are illustrative (inside the range fitted to lipid agents, not one agent).
% Runtime : about 15 s in Octave (RK4 in an interpreted loop), a few s in MATLAB.
clear; close all; rng(1);

%% Liquid and gas constants (water at 20 C)
rho = 1000;        % kg/m^3
c   = 1480;        % m/s
mu  = 1e-3;        % Pa s
P0  = 101325;      % Pa
sw  = 0.072;       % N/m, water surface tension
rObs = 0.01;       % m, observation distance for the emitted pressure

%% Cases: {label, R0 (m), shell (0/1), gamma, s0 (N/m), f (Hz or 0 = 2 f0), P (Pa)}
chi = 0.5;  ks = 5e-9;            % shell elasticity (N/m) and dilatational viscosity (kg/s)
cases = { 'free, f = 2 MHz,  50 kPa',  2e-6, 0, 1.40, 0,    2e6,  50e3;
          'free, f = 2 MHz, 200 kPa',  2e-6, 0, 1.40, 0,    2e6, 200e3;
          'free, f = 2 f0,  100 kPa',  2e-6, 0, 1.40, 0,    0,   100e3;
          'free, f = 2 f0,  200 kPa',  2e-6, 0, 1.40, 0,    0,   200e3;
          'shell elastic s0=0.02, 3 MHz, 100 kPa', 2e-6, 1, 1.07, 0.02, 3e6, 100e3;
          'shell buckled s0=0, 3 MHz, 50 kPa',     2e-6, 1, 1.07, 0,    3e6,  50e3 };
ncyc = 30;        % cycles in the tone burst (sin^2 ramps of one cycle)
nper = 200;       % RK4 steps per acoustic period (split further when R moves fast)

nc = size(cases, 1);
res = zeros(nc, 7);
keep = cell(nc, 1);
fprintf('%-40s %7s %6s %6s %7s %7s %7s %8s\n', 'case', 'f0 MHz', 'Rmax', 'Rmin', 'H2 dB', 'H3 dB', 'f/2 dB', 'ps,max Pa');
for ic = 1:nc
  R0 = cases{ic, 2}; shell = cases{ic, 3}; g = cases{ic, 4}; s0 = cases{ic, 5};
  %% Shell law and equilibrium gas pressure
  if shell
    Rb = R0 / sqrt(1 + s0 / chi);  Rr = Rb * sqrt(1 + sw / chi);
    sig = @(R) (R >= Rb & R <= Rr) .* chi .* (R.^2 / Rb^2 - 1) + (R > Rr) .* sw;
    sig0 = s0; kss = ks;
    K = 3*g*P0 + 2*(3*g + 1)*sig0/R0 + 4*chi/R0;      % linear stiffness (elastic state)
  else
    sig = @(R) sw + 0*R; sig0 = sw; kss = 0;
    K = 3*g*P0 + 2*(3*g - 1)*sw/R0;
  end
  f0 = sqrt(K / rho) / (2*pi*R0);
  f = cases{ic, 6}; if f == 0, f = 2 * f0; end
  P = cases{ic, 7};
  Pg0 = P0 + 2 * sig0 / R0;
  T = ncyc / f; tr = 1 / f;
  env = @(t) (t > 0 & t < T) .* sin(pi/2 * min(1, min(t / tr, (T - t) / tr))).^2;
  pac = @(t) P * env(t) .* sin(2*pi*f*t);
  acc = @(t, R, U) ((Pg0 * (R0 ./ R).^(3*g) .* (1 - 3*g*U/c) - P0 - 2*sig(R)./R ...
                    - 4*mu*U./R - 4*kss*U./R.^2 - pac(t)) / rho - 1.5*U.^2) ./ R;
  %% RK4 integration
  dt = 1 / (f * nper); n = round((T + 4/f) / dt);
  t = (0:n) * dt; R = zeros(1, n + 1); U = zeros(1, n + 1); A = zeros(1, n + 1);
  R(1) = R0; A(1) = acc(0, R0, 0); Rk = R0; Uk = 0; tk = 0;
  for i = 1:n
    m = max(1, min(64, ceil(abs(Uk) * dt / Rk / 0.005)));   % sub-steps when R moves > 0.5 % per step
    h = dt / m;
    for s = 1:m
      k1R = Uk;              k1U = acc(tk, Rk, Uk);
      k2R = Uk + h/2*k1U;    k2U = acc(tk + h/2, Rk + h/2*k1R, Uk + h/2*k1U);
      k3R = Uk + h/2*k2U;    k3U = acc(tk + h/2, Rk + h/2*k2R, Uk + h/2*k2U);
      k4R = Uk + h*k3U;      k4U = acc(tk + h,   Rk + h*k3R,   Uk + h*k3U);
      Rk = Rk + h/6*(k1R + 2*k2R + 2*k3R + k4R);
      Uk = Uk + h/6*(k1U + 2*k2U + 2*k3U + k4U);
      tk = tk + h;
    end
    R(i+1) = Rk; U(i+1) = Uk; A(i+1) = acc(tk, Rk, Uk);
  end
  ps = rho * (R.^2 .* A + 2 * R .* U.^2) / rObs;     % emitted pressure at rObs (Pa)
  %% Spectrum of the steady part (30 % to 95 % of the burst), Hann window, 8x zero padding
  idx = find(t > 0.3*T & t < 0.95*T); x = ps(idx) - mean(ps(idx)); L = numel(x);
  w = 0.5 - 0.5*cos(2*pi*(0:L-1)/(L-1));
  nf = 2^nextpow2(8*L); X = abs(fft(x .* w, nf)); fa = (0:nf-1) / (nf*dt);
  lev = @(ff) max(X(max(1, round(ff/fa(2)) - 2):round(ff/fa(2)) + 4));
  L1 = lev(f);
  res(ic, :) = [f0/1e6, max(R)/R0, min(R)/R0, 20*log10(lev(2*f)/L1), 20*log10(lev(3*f)/L1), 20*log10(lev(f/2)/L1), max(ps)];
  fprintf('%-40s %7.2f %6.2f %6.2f %7.1f %7.1f %7.1f %8.1f\n', cases{ic, 1}, res(ic, :));
  keep{ic} = struct('t', t, 'R', R / R0, 'ps', ps, 'f', f, 'fa', fa(1:nf/2), 'X', 20*log10(X(1:nf/2)/L1 + 1e-12));
end

%% Compression-only index for the buckled shell: (R0 - Rmin)/(Rmax - R0)
fprintf('\nCompression/expansion ratio (R0 - Rmin)/(Rmax - R0): elastic %.2f, buckled %.2f\n', ...
  (1 - res(5, 3)) / (res(5, 2) - 1), (1 - res(6, 3)) / (res(6, 2) - 1));

%% Step-size check on the most violent case (free, 2 MHz, 200 kPa): already converged at 200/period
fprintf('Step check: Rmin/R0 = %.4f at %d steps/period (compare 0.415 at 400 and 800 steps/period)\n', res(2, 3), nper);

%% Figures
figure('Position', [50 50 900 640]);
sel = [1 2 4 6];
for k = 1:4
  s = keep{sel(k)};
  subplot(4, 2, 2*k - 1); plot(s.t * 1e6, s.R, 'b'); grid on;
  ylabel('R / R_0'); title(cases{sel(k), 1});
  if k == 4, xlabel('time (\mus)'); end
  subplot(4, 2, 2*k); plot(s.fa / s.f, s.X, 'r'); grid on; xlim([0 4]); ylim([-80 5]);
  ylabel('dB re f'); title('emitted spectrum');
  if k == 4, xlabel('frequency / drive frequency'); end
end
