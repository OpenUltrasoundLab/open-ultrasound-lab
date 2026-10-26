%% m28_resonance_shell.m  --  Minnaert, surface tension, shell stiffness, damping, cross-section
% Purpose : (1) evaluate the Minnaert constant f0*R0 for the usual assumptions;
%           (2) add surface tension and a Marmottant shell to the linear resonance;
%           (3) verify the linear formula by a small-amplitude RK4 ring-down;
%           (4) compute damping (radiation, liquid viscosity, shell viscosity) and the
%               scattering cross-section, and compare with a rigid sphere.
% Linear  : rho w0^2 R0^2 = 3 g P0 + 2(3g - 1) sigma/R0                 (free bubble)
%           rho w0^2 R0^2 = 3 g P0 + 2(3g + 1) s0/R0 + 4 chi/R0         (shell, elastic state)
%           sigma_s = 4 pi R0^2 / [((f0/f)^2 - 1)^2 + delta^2],
%           delta = w R0/c + 4 mu/(rho w R0^2) + 4 ks/(rho w R0^3)
% Change  : chi, ks, s0, gamma, the radius list.
% Expect  : f0*R0 = 3.29 MHz um (air, 1 atm, adiabatic), 3.26 at 1 bar, 2.78 isothermal;
%           a free 1 um air bubble resonates at 4.74 MHz (not 3.29) because of surface
%           tension; a 2 um shelled bubble (g = 1.07) at 2.99 MHz versus 1.44 MHz
%           (Minnaert, same gas) and 1.75 MHz (free, with tension); ring-down
%           frequencies agree with sqrt(w0^2 - (b/2)^2) to better than 0.1 %.
% Provenance: closed-form linear formulas checked by an RK4 ring-down. Shell values chi = 0.5 N/m, s0 = 0.02 N/m, ks = 5e-9 kg/s,
%           g = 1.07 are illustrative (inside the range fitted to lipid agents, not one agent).
% Runtime : about 2 s.
clear; close all; rng(1);

rho = 1000; c = 1480; mu = 1e-3; P0 = 101325; sw = 0.072;
chi = 0.5; ks = 5e-9; s0 = 0.02;          % shell: N/m, kg/s, N/m

%% 1. Minnaert constant (f0 * R0 in m/s = MHz * um)
mc = @(g, p, rr) sqrt(3 * g * p / rr) / (2*pi);
fprintf('Minnaert f0*R0 (MHz um): air 1 atm adiabatic %.3f | 1 bar, rho 1000: %.3f | isothermal %.3f | g = 1.07 (C3F8/C4F10-like) %.3f | SF6 g = 1.095: %.3f\n', ...
  mc(1.4, P0, 998), mc(1.4, 1e5, 1000), mc(1.0, P0, rho), mc(1.07, P0, rho), mc(1.095, P0, rho));

%% 2. Resonance versus radius for four models
R0s = [0.5 1 1.5 2 3 4 5] * 1e-6;
fM  = @(R, g) sqrt(3*g*P0/rho) ./ (2*pi*R);
fT  = @(R, g) sqrt((3*g*P0 + 2*(3*g - 1)*sw./R) / rho) ./ (2*pi*R);
fS  = @(R, g) sqrt((3*g*P0 + 2*(3*g + 1)*s0./R + 4*chi./R) / rho) ./ (2*pi*R);
fprintf('\n R0 (um) | Minnaert air | + tension air | free g=1.07 + tension | shell g=1.07 (MHz)\n');
for R = R0s
  fprintf('  %4.1f   |   %6.2f     |    %6.2f     |        %6.2f         |   %6.2f\n', R*1e6, fM(R, 1.4)/1e6, fT(R, 1.4)/1e6, fT(R, 1.07)/1e6, fS(R, 1.07)/1e6);
end

%% 3. Ring-down check with RK4 (release from R = 1.002 R0, no drive)
fprintf('\nRing-down check (damped frequency, MHz)\n');
cfg = {'free air', 1.4, 0; 'shell', 1.07, 1};
for ic = 1:2
  g = cfg{ic, 2}; shell = cfg{ic, 3};
  for R0 = [1 2 3] * 1e-6
    if shell
      Rb = R0 / sqrt(1 + s0/chi); Rr = Rb * sqrt(1 + sw/chi);
      sig = @(R) (R >= Rb & R <= Rr) .* chi .* (R.^2/Rb^2 - 1) + (R > Rr) .* sw; sig0 = s0; kk = ks;
      w0 = 2*pi*fS(R0, g);
    else
      sig = @(R) sw + 0*R; sig0 = sw; kk = 0; w0 = 2*pi*fT(R0, g);
    end
    Pg0 = P0 + 2*sig0/R0;
    acc = @(R, U) ((Pg0*(R0./R).^(3*g).*(1 - 3*g*U/c) - P0 - 2*sig(R)./R - 4*mu*U./R - 4*kk*U./R.^2)/rho - 1.5*U.^2)./R;
    b = (4*mu/R0 + 4*kk/R0^2 + 3*g*Pg0/c) / (rho*R0);      % linear damping coefficient (1/s)
    wd = sqrt(w0^2 - (b/2)^2);
    dt = 2*pi/w0/400; n = 2400; x = zeros(1, n); Rk = 1.002*R0; Uk = 0;
    for i = 1:n
      k1R = Uk; k1U = acc(Rk, Uk);
      k2R = Uk + dt/2*k1U; k2U = acc(Rk + dt/2*k1R, Uk + dt/2*k1U);
      k3R = Uk + dt/2*k2U; k3U = acc(Rk + dt/2*k2R, Uk + dt/2*k2U);
      k4R = Uk + dt*k3U;   k4U = acc(Rk + dt*k3R,   Uk + dt*k3U);
      Rk = Rk + dt/6*(k1R + 2*k2R + 2*k3R + k4R); Uk = Uk + dt/6*(k1U + 2*k2U + 2*k3U + k4U);
      x(i) = Rk - R0;
    end
    zc = find(x(1:end-1) < 0 & x(2:end) >= 0);            % upward zero crossings
    zf = zc - x(zc) ./ (x(zc+1) - x(zc));                   % sub-sample interpolation
    fmeas = 1 / (mean(diff(zf(1:min(4, end)))) * dt);
    fprintf('  %-8s R0 = %.0f um: formula f0 = %.3f, damped sqrt(w0^2-(b/2)^2) = %.3f, RK4 = %.3f, Q = %.1f\n', ...
      cfg{ic, 1}, R0*1e6, w0/2/pi/1e6, wd/2/pi/1e6, fmeas/1e6, w0/b);
  end
end

%% 4. Damping and cross-section at resonance (physical radiation damping w R0 / c)
fprintf('\nDamping at resonance and cross-section (free bubble g = 1.4 with tension; shell g = 1.07)\n');
fprintf(' R0 (um) | model | f0 MHz | d_rad  d_vis  d_shell | Q   | sigma_s/(pi R0^2) at f0\n');
for R0 = [1 2 3] * 1e-6
  for shell = 0:1
    if shell, f0 = fS(R0, 1.07); kk = ks; lab = 'shell'; else, f0 = fT(R0, 1.4); kk = 0; lab = 'free '; end
    w = 2*pi*f0; dr = w*R0/c; dv = 4*mu/(rho*w*R0^2); dsh = 4*kk/(rho*w*R0^3); d = dr + dv + dsh;
    fprintf('  %3.0f    | %s | %5.2f  | %.3f  %.3f  %.3f  | %4.1f | %6.0f\n', R0*1e6, lab, f0/1e6, dr, dv, dsh, 1/d, 4/d^2);
  end
end
kappaRatio = rho*c^2/(1.4*P0);
fprintf('\nGas/water compressibility ratio (air, adiabatic): %.0f; bubble / rigid sphere cross-section far below resonance: (4/7) ratio^2 = %.2e\n', ...
  kappaRatio, 4/7*kappaRatio^2);

%% 5. Figures: resonance versus radius, and cross-section versus frequency for R0 = 2 um
Rg = linspace(0.5, 5, 200) * 1e-6;
figure('Position', [50 50 900 380]);
subplot(1, 2, 1);
loglog(Rg*1e6, fM(Rg, 1.4)/1e6, 'k--', Rg*1e6, fT(Rg, 1.4)/1e6, 'b', Rg*1e6, fT(Rg, 1.07)/1e6, 'c', Rg*1e6, fS(Rg, 1.07)/1e6, 'r'); grid on;
xlabel('R_0 (\mum)'); ylabel('f_0 (MHz)'); title('Linear resonance'); xlim([0.5 5]);
legend('Minnaert, air', 'air + tension', '\gamma = 1.07 + tension', 'lipid shell', 'Location', 'southwest');
subplot(1, 2, 2);
fg = linspace(0.2, 10, 600) * 1e6; R0 = 2e-6;
for shell = 0:1
  if shell, f0 = fS(R0, 1.07); kk = ks; st = 'r'; else, f0 = fT(R0, 1.4); kk = 0; st = 'b'; end
  w = 2*pi*fg; d = w*R0/c + 4*mu./(rho*w*R0^2) + 4*kk./(rho*w*R0^3);
  ss = 4*pi*R0^2 ./ (((f0./fg).^2 - 1).^2 + d.^2);
  semilogy(fg/1e6, ss/(pi*R0^2), st); hold on;
end
k = 2*pi*fg/c; semilogy(fg/1e6, 7/9*(k*R0).^4, 'k--');
grid on; xlabel('f (MHz)'); ylabel('\sigma_s / \pi R_0^2'); title('Scattering cross-section, R_0 = 2 \mum');
legend('free (air)', 'lipid shell', 'rigid sphere', 'Location', 'southeast'); ylim([1e-10 1e4]);
