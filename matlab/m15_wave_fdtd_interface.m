%% m15_wave_fdtd_interface.m
% Purpose : 1D linear acoustic FDTD (staggered pressure / particle velocity, leapfrog)
%           launched at an impedance step and through a quarter-wave layer.
%           Measured reflected / transmitted ENERGY fractions are compared with
%           R = (Z2-Z1)/(Z2+Z1) and 1 - R^2 (band-averaged for the layer case).
%           A second cell computes oblique-incidence energy partition at a
%           soft-tissue / skull-bone interface, including mode conversion to shear.
% Change  : the cases table (c, rho of each medium, layer), f0, FBW, PPW.
% Expect  : fat->muscle  R^2 = 0.0086 (R = 0.093, -20.7 dB)
%           tissue->skull R^2 = 0.288 ; tissue->air-like R^2 = 0.997
%           PZT (30 MRayl)->tissue R^2 = 0.81 bare, 0.18 band-averaged with a lambda/4 layer of
%           sqrt(Z1 Z2) (zero at f0 only; 60 % pulse)
%           skull: critical angle 33.1 deg, > 90 % of the energy enters as shear at 40-70 deg
% Runtime : about 3 s in Octave.
clear; close all; rng(1);

%% Parameters
f0  = 3e6;          % centre frequency (Hz)
FBW = 0.6;          % -6 dB fractional bandwidth of the Gaussian pulse
PPW = 24;           % points per wavelength at f0 in the slowest medium
sf  = FBW*f0/(2*sqrt(2*log(2)));   % spectral std (Hz)
st  = 1/(2*pi*sf);                 % envelope std (s)
pulse = @(t) exp(-t.^2/(2*st^2)).*sin(2*pi*f0*t);

% name, c1, rho1, c2, rho2, layer thickness (m), c_layer, rho_layer
cases = { 'fat -> muscle',            1450,  950, 1580, 1050, 0, 0, 0;
          'soft tissue -> skull',     1540, 1050, 2820, 1900, 0, 0, 0;
          'soft tissue -> air-like',  1540, 1050, 1000,  1.2, 0, 0, 0;
          'PZT -> tissue, bare',      4000, 7500, 1540, 1050, 0, 0, 0;
          'PZT -> tissue, lambda/4',  4000, 7500, 1540, 1050, 2500/(4*f0), 2500, NaN };
Z1Z2 = 4000*7500*1540*1050; cases{5,8} = sqrt(Z1Z2)/2500;   % layer impedance sqrt(Z1 Z2)

fprintf('%-26s %8s %8s %8s %8s %9s %9s\n', 'case', 'R th', 'R meas', 'T th', 'T meas', 'Er th', 'Er meas');
snap = {};
for ic = 1:size(cases, 1)
  c1 = cases{ic,2}; r1 = cases{ic,3}; c2 = cases{ic,4}; r2 = cases{ic,5};
  dL = cases{ic,6}; cL = cases{ic,7}; rL = cases{ic,8};
  cmin = min([c1 c2]); cmax = max([c1 c2]); if dL > 0, cmin = min(cmin, cL); cmax = max(cmax, cL); end
  dx = cmin/f0/PPW; if dL > 0, dx = dL/ceil(dL/dx); end   % layer = whole number of cells
  dt = 0.95*dx/cmax;
  xI = (round(15e-3/dx) - 0.5)*dx;     % interface position (m), midway between p nodes
  x0 = xI - 9e-3;                      % initial pulse centre
  xA = xI - 4e-3; xB = xI + dL + 2e-3; % receivers in medium 1 and medium 2
  tEnd = (xA - x0 + 8e-3)/c1 + 3e-6;   % reflected pulse (and layer ring-down) past receiver A
  xEnd = xB + c2*tEnd + 2e-3;
  x  = (0:dx:xEnd)';  nx = numel(x);  xh = x + dx/2;      % p at x, u at x + dx/2
  c = c2*ones(nx,1); rho = r2*ones(nx,1);
  c(x < xI) = c1; rho(x < xI) = r1;
  if dL > 0, in = (x >= xI) & (x < xI + dL); c(in) = cL; rho(in) = rL; end
  K = rho.*c.^2;                                  % bulk modulus at p nodes
  rh = rho; rh(1:end-1) = 0.5*(rho(1:end-1) + rho(2:end));   % density at u nodes
  % initial right-going pulse: p(x,0) = s(-(x-x0)/c1); u = p/Z1 half a time step later
  p = pulse(-(x - x0)/c1);                 p(x >= xI) = 0;
  u = pulse(-(xh - x0)/c1 + dt/2)/(r1*c1); u(xh >= xI) = 0;   % u lives at t = +dt/2 (leapfrog)
  nt = ceil(tEnd/dt); iA = round(xA/dx) + 1; iB = round(xB/dx) + 1;
  pA = zeros(nt,1); pB = zeros(nt,1);
  for n = 1:nt
    p(2:end) = p(2:end) - dt*K(2:end).*(u(2:end) - u(1:end-1))/dx;          % dp/dt = -K du/dx
    u(1:end-1) = u(1:end-1) - dt./rh(1:end-1).*(p(2:end) - p(1:end-1))/dx;  % du/dt = -(1/rho) dp/dx
    pA(n) = p(iA); pB(n) = p(iB);
    if n == round((xI - x0 + 3e-3)/c1/dt), snap{ic} = [x p]; end %#ok<SAGROW>
  end
  t = (1:nt)'*dt; tm = (xA - x0 + 4e-3)/c1;          % split incident / reflected at receiver A
  inc = pA.*(t < tm); ref = pA.*(t >= tm);
  % pressure coefficients by projection onto the incident waveform (shifted to each arrival)
  shiftTo = @(sig, d) interp1(t, sig, t - d, 'linear', 0);
  incR = shiftTo(inc, 8e-3/c1);                              % incident shape at the reflected arrival
  incT = shiftTo(inc, (xI - xA)/c1 + dL*(dL > 0)/max(cL,1) + (xB - xI - dL)/c2);
  % allow a small residual shift (the material jump sits between grid nodes): best-aligned projection
  proj = @(sig, ref0) max(arrayfun(@(d) abs(sum(sig.*shiftTo(ref0, d))), (-40:40)*dt/4)) / sum(ref0.^2);
  Rm = sign(sum(ref.*incR))*proj(ref, incR); Tm = proj(pB, incT);
  Erm = sum(ref.^2)/sum(inc.^2);                             % reflected energy flux ratio (same medium)
  Rth = (r2*c2 - r1*c1)/(r2*c2 + r1*c1); Tth = 1 + Rth;
  if dL == 0
    Erth = Rth^2;
  else  % band-averaged energy reflection of the layer, weighted by the pulse power spectrum
    f = linspace(f0 - 4*sf, f0 + 4*sf, 801); f = f(f > 0); w = exp(-(f - f0).^2/sf^2);
    Za = r1*c1; Zb = r2*c2; Zl = rL*cL; kd = 2*pi*f*dL/cL;
    Zin = Zl*(Zb + 1i*Zl*tan(kd))./(Zl + 1i*Zb*tan(kd));    % input impedance seen from medium 1
    Rf = (Zin - Za)./(Zin + Za); Erth = sum(w.*abs(Rf).^2)/sum(w);
    Rth = NaN; Tth = NaN; Rm = NaN; Tm = NaN;                % waveform changes: compare energy only
  end
  fprintf('%-26s %8.4f %8.4f %8.4f %8.4f %9.4f %9.4f\n', cases{ic,1}, Rth, Rm, Tth, Tm, Erth, Erm);
end
fprintf('(T = 1 + R is a pressure ratio and can exceed 1; the energy fractions are R^2 and 1 - R^2.)\n');

%% Figure: snapshots
figure('Position', [50 50 900 640]);
for ic = 1:size(cases, 1)
  subplot(size(cases,1), 1, ic); s = snap{ic};
  plot(s(:,1)*1e3, s(:,2), 'b'); hold on; yl = [-1.6 1.6]; ylim(yl);
  plot([15 15], yl, 'k--'); if cases{ic,6} > 0, plot([15 15]+cases{ic,6}*1e3, yl, 'k:'); end
  xlim([0 35]); ylabel('p / p_{inc}'); title(cases{ic,1});
end
xlabel('x (mm)');

%% Oblique incidence: fluid (tissue) -> isotropic solid (skull), energy partition
% Potentials: fluid phi1 = e^{i(kx x + kz1 z)} + R e^{i(kx x - kz1 z)};
% solid phi = A e^{i(kx x + kzL z)}, psi = B e^{i(kx x + kzS z)}. Boundary conditions at z = 0:
% u_z continuous, sigma_zz = -p, sigma_xz = 0.
r1 = 1050; c1 = 1540; r2 = 1900; cLb = 2820; cSb = 1500;   % White, Clement & Hynynen 2006 speeds (1 MHz)
lam = r2*(cLb^2 - 2*cSb^2); mu = r2*cSb^2; w = 1;
th = (0:1:89)*pi/180; ER = zeros(size(th)); EL = ER; ES = ER;
for i = 1:numel(th)
  k1 = w/c1; kL = w/cLb; kS = w/cSb; kx = k1*sin(th(i));
  kz1 = sqrt(k1^2 - kx^2); kzL = sqrt(kL^2 - kx^2 + 0i); kzS = sqrt(kS^2 - kx^2 + 0i);
  if imag(kzL) < 0, kzL = -kzL; end
  if imag(kzS) < 0, kzS = -kzS; end
  M = [1i*kz1, 1i*kzL, 1i*kx;
       r1*w^2, -lam*kL^2 - 2*mu*kzL^2, -2*mu*kx*kzS;
       0, -2*kx*kzL, kzS^2 - kx^2];
  v = M \ [1i*kz1; -r1*w^2; 0];
  ER(i) = abs(v(1))^2;
  EL(i) = r2*cLb*kL^2*real(kzL/kL)/(r1*c1*k1^2*cos(th(i)))*abs(v(2))^2;
  ES(i) = r2*cSb*kS^2*real(kzS/kS)/(r1*c1*k1^2*cos(th(i)))*abs(v(3))^2;
end
fprintf('\nTissue -> skull (cL %d, cS %d m/s): longitudinal critical angle %.1f deg\n', cLb, cSb, asind(c1/cLb));
fprintf(' angle  E_refl  E_long  E_shear   sum\n');
for a = [0 20 30 40 50 60 70 80]
  i = a + 1; fprintf('%5d  %6.3f  %6.3f  %6.3f  %6.4f\n', a, ER(i), EL(i), ES(i), ER(i) + EL(i) + ES(i));
end
figure('Position', [80 80 640 360]);
plot(th*180/pi, ER, 'k', th*180/pi, EL, 'b', th*180/pi, ES, 'r', 'LineWidth', 1.4);
xlabel('incidence angle (deg)'); ylabel('energy fraction'); grid on;
legend('reflected', 'longitudinal in bone', 'shear in bone', 'Location', 'east');
title('Soft tissue -> skull: mode conversion');
