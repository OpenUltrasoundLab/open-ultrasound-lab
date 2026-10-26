%% m18_cmut_vs_pzt_bandwidth.m
% Pulse-echo bandwidth of a PZT-5H stack (Mason model, backing, 0-2 quarter-wave
% matching layers) versus a CMUT element (lumped piston membrane, water load).
% APPROXIMATE MODEL: 1D, plane-wave radiation resistance rho*c*A for both,
% no diffraction, no lens, lossless layers; the CMUT membrane is one piston mode
% (no higher membrane modes, no crosstalk), with an optional cell-level water
% added mass (8/3) rho a^3 per cell.
%
% What it shows: the PZT stack is a resonator that matching layers broaden;
%   the CMUT is an over-damped mass-spring whose band runs from k_eff/(rho c)
%   to rho c / m, so its upper edge is set by the moving mass (membrane + water).
% What to change: ZB (backing), nL (layers), bias fraction Vf, a, h.
% Expected output (defaults): PZT 2 layers, ZB = 5 MRayl: FBW 63 %;
%   CMUT a = 20 um, h = 1 um, 0.9 Vpi with added mass: FBW 178 % (model upper
%   bound; measured immersion CMUT arrays report about 90-105 %).
% Runtime: about 2 s.
clear; close all; rng(1);
e0 = 8.854e-12; rw = 1000; cw = 1500;
A = 0.27e-3*5e-3; Rs = 50; ZL = rw*cw*A;     % element 0.27 mm x 5 mm, 50 ohm source and receiver
f = (1:2000)'*12.5e3; w = 2*pi*f;              % 12.5 kHz to 25 MHz

%% PZT-5H stack, Mason model (Ma et al. 2015: v 4580 m/s, kt 0.51, epsS 1470; rho 7500)
v = 4580; rho = 7500; kt = 0.51; er = 1470; f0 = 5e6;
d = v/(2*f0); Z0 = rho*v*A; C0 = er*e0*A/d; hh = kt*sqrt(rho*v^2/(er*e0));
Zp = rho*v/1e6; cm = 2500;
SP = zeros(numel(f), 3);
for nL = 0:2
  if nL == 2, Zm = [Zp^(4/7)*1.5^(3/7), Zp^(1/7)*1.5^(6/7)];   % Desilets (Module 17)
  elseif nL == 1, Zm = (Zp*1.5^2)^(1/3); else, Zm = []; end
  ZB = 5e6*A;
  for k = 1:numel(f)
    th = w(k)*d/v; aa = Z0/(1i*tan(th)); bb = Z0/(1i*sin(th)); cc = hh/(1i*w(k));
    Zin = ZL; T = eye(2);
    for m = numel(Zm):-1:1                      % from the water side towards the plate
      Zc = Zm(m)*1e6*A; t = w(k)/(4*f0);         % quarter wave at f0 (t = omega d / c, d = c/(4 f0))
      Zin = Zc*(Zin + 1i*Zc*tan(t))/(Zc + 1i*Zin*tan(t));
      T = [cos(t), 1i*Zc*sin(t); 1i*sin(t)/Zc, cos(t)] * T;   % plate-side layer first after the loop
    end
    Mx = [aa+ZB, bb, cc; bb, aa+Zin, cc; cc, cc, 1/(1i*w(k)*C0)+Rs];
    x = Mx\[0;0;1];
    Fv = [-Zin*x(2); -x(2)];                     % force and outward velocity at the plate front
    out = T\Fv;                                  % at the water
    SP(k, nL+1) = 2*sqrt(Rs/ZL)*out(1);
  end
end

%% CMUT element: lumped piston (average-displacement spring, 1.8 rho h modal mass)
E = 169e9; nus = 0.22; rs = 2330; a = 20e-6; hm = 1e-6; g0 = 0.15e-6; ti = 0.1e-6; eri = 7.5; eta = 0.8;
Dm = E*hm^3/(12*(1-nus^2)); ge = g0 + ti/eri;
k = 192*Dm/a^4*A*eta; mm = 1.8*rs*hm*A*eta; madd = (8/3)*rw*a^3/(pi*a^2)*A*eta;
Vpi = sqrt(8*k*ge^3/(27*e0*eta*A));
xi = linspace(0, 1/3, 20001); Vr = sqrt(27/4*xi.*(1-xi).^2);
SC = zeros(numel(f), 3); Vfs = [0.7 0.9 0.9]; adds = [1 1 0];
for j = 1:3
  Vf = Vfs(j); x = interp1(Vr, xi, Vf)*ge; Vdc = Vf*Vpi;
  C0c = e0*eta*A/(ge - x); n = C0c*Vdc/(ge - x); keff = k*(1 - 2*(x/ge)/(1 - x/ge));
  mt = mm + adds(j)*madd;
  Zmech = 1i*w*mt + keff./(1i*w);
  Aa = 1/n; B = Zmech/n; Cc = 1i*w*C0c/n; Dd = 1i*w.*C0c.*Zmech/n + n;
  SC(:, j) = 2*sqrt(Rs*ZL)./(Aa*ZL + B + Rs*(Cc*ZL + Dd));
end
fprintf('CMUT: f0(vac) = %.2f MHz, Vpi = %.1f V, membrane mass %.2g, cell added mass %.2g kg/m^2\n', ...
  10.2158/(2*pi*a^2)*sqrt(Dm/(rs*hm))/1e6, Vpi, 1.8*rs*hm, (8/3)*rw*a/pi);

%% Bandwidth metrics on the pulse-echo response H = S21^2
labels = {'PZT, no matching', 'PZT, 1 layer', 'PZT, 2 layers', 'CMUT 0.7 Vpi', 'CMUT 0.9 Vpi', 'CMUT 0.9 Vpi, no added mass'};
H = [SP.^2, SC.^2];
for j = 1:size(H,2)
  m = abs(H(:,j)); [pk, ip] = max(m);
  il = find(m(1:ip) < pk/2, 1, 'last'); ih = ip - 1 + find(m(ip:end) < pk/2, 1, 'first');
  if isempty(il), fl = f(1); else, fl = f(il) + (pk/2 - m(il))/(m(il+1) - m(il))*(f(il+1) - f(il)); end     % linear interpolation of the -6 dB edges
  if isempty(ih), fh = f(end); else, fh = f(ih-1) + (m(ih-1) - pk/2)/(m(ih-1) - m(ih))*(f(ih) - f(ih-1)); end
  fprintf('%-28s fc = %5.2f MHz, -6 dB band %5.2f-%5.2f MHz, FBW = %4.0f %%, peak %6.1f dB\n', labels{j}, ...
    (fl+fh)/2/1e6, fl/1e6, fh/1e6, 100*(fh-fl)/((fh+fl)/2), 20*log10(pk));
end

%% Impulse responses: inverse FFT of H, receive low-pass (cos^2 roll-off 9-12 MHz, as in Module 17), 0.5 us delay
Nf = 4096; fs = Nf*12.5e3; Hfull = zeros(Nf, 2);
Hfull(2:numel(f)+1, 1) = H(:,3); Hfull(2:numel(f)+1, 2) = H(:,5);
win = (f < 9e6) + (f >= 9e6 & f < 12e6).*cos(pi/2*(f - 9e6)/3e6).^2;
win = win.*exp(-1i*w*0.5e-6);
Hfull(2:numel(f)+1, :) = Hfull(2:numel(f)+1, :).*[win win];
Hfull(Nf:-1:Nf/2+2, :) = conj(Hfull(2:Nf/2, :));
hpe = real(ifft(Hfull));
t = (0:Nf-1)/fs*1e6;

figure('Position', [80 80 1000 380]);
subplot(1,2,1);
plot(f/1e6, 20*log10(abs(H(:,3))/max(abs(H(:,3)))), 'b-', f/1e6, 20*log10(abs(H(:,5))/max(abs(H(:,5)))), 'r-', 'LineWidth', 2);
grid on; ylim([-30 2]); xlabel('frequency (MHz)'); ylabel('pulse-echo |H| (dB, normalized)');
legend('PZT-5H, 2 layers, Z_B 5 MRayl', 'CMUT 0.9 V_{PI}', 'Location', 'south'); title('Approximate 1D model');
subplot(1,2,2);
plot(t, hpe(:,1)/max(abs(hpe(:,1))), 'b-', t, hpe(:,2)/max(abs(hpe(:,2))) - 2.2, 'r-');
xlim([0.3 2.0]); xlabel('time (\mus)'); ylabel('normalized'); title('Impulse responses, 9-12 MHz receive low-pass'); grid on;
