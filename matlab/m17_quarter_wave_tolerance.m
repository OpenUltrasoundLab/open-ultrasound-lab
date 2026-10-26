%% m17_quarter_wave_tolerance.m  --  one quarter-wave layer: bandwidth, tolerance, loss, bondline
% Purpose : semi-infinite piezo (Zp) radiating into tissue (ZL) through one
%           layer. Checks Z1 = sqrt(Zp ZL), the -3 dB transmission band against
%           the closed-form quarter-wave-transformer formula, and what thickness
%           error, layer attenuation and a thin glue bondline do.
% Change  : Zp (a 1-3 composite is 10-25 MRayl), ZL, bond thickness, alpha.
% Expect  : T(no layer) = 0.161 (-7.9 dB); -3 dB band 0.711-1.289 f0 = 57.7 %
%           (formula 57.7 %); +-10 % thickness error: T(f0) = 0.887, band
%           centre moves by -+10 %; alpha = 0.5 dB/mm/MHz absorbs 0.73 dB one
%           way at ANY f0 (alpha ~ f, d ~ 1/f); 5 um bond (3 MRayl): -0.5 dB at
%           5 MHz, -4.7 dB at 20 MHz. Runtime ~ 1 s.
clear; close all; rng(1);

%% Media (MRayl) and helpers
Zp = 34.2; ZL = 1.5; Z1 = sqrt(Zp*ZL);
f = linspace(0.001, 2, 4001);                     % frequency / f0
% ABCD of a line with electrical angle th (rad), impedance Z, loss factor a*d (Np)
lineM = @(Z, th, ad) [cosh(ad + 1i*th), Z*sinh(ad + 1i*th); sinh(ad + 1i*th)/Z, cosh(ad + 1i*th)];
Tpow = @(Zin) 1 - abs((Zin - Zp)./(Zin + Zp)).^2;  % power transmission from Zp into the stack

%% 1. No layer and the ideal quarter-wave layer
T0 = 4*Zp*ZL/(Zp + ZL)^2;
th = pi/2*f;                                      % electrical length, pi/2 at f0
Zin = Z1*(ZL + 1i*Z1*tan(th))./(Z1 + 1i*ZL*tan(th));
T1 = Tpow(Zin);
band = f(T1 >= 0.5);
Gm = sqrt(0.5);                                   % |Gamma| at the -3 dB transmission edge
bwFormula = 2 - 4/pi*acos(2*Gm*sqrt(Zp*ZL)/(abs(ZL - Zp)*sqrt(1 - Gm^2)));
fprintf('No layer: T = %.3f (%.1f dB), R_p = %.3f\n', T0, 10*log10(T0), (ZL - Zp)/(ZL + Zp));
fprintf('Z1 = sqrt(Zp ZL) = %.2f MRayl; Zin(f0) = Z1^2/ZL = %.2f MRayl\n', Z1, Z1^2/ZL);
fprintf('-3 dB transmission band: %.3f to %.3f f0 = %.1f %% (closed form %.1f %%)\n', band(1), band(end), 100*(band(end) - band(1)), 100*bwFormula);
r1 = (Z1 - Zp)/(Z1 + Zp); r2 = (ZL - Z1)/(ZL + Z1);
fprintf('Two interface reflections at Z1 = sqrt: r1 = %.3f, r2 = %.3f (equal, cancel at f0)\n', r1, r2);

%% 2. Thickness error: scale the electrical length by s
sv = [0.9 0.95 1.05 1.1]; Ts = zeros(numel(sv), numel(f));
for k = 1:numel(sv)
  tk = pi/2*f*sv(k); Ts(k,:) = Tpow(Z1*(ZL + 1i*Z1*tan(tk))./(Z1 + 1i*ZL*tan(tk)));
  [~, ip] = max(Ts(k,:));
  [~, i0] = min(abs(f - 1));
  fprintf('thickness x %.2f: T(f0) = %.3f (%.2f dB), peak at %.3f f0\n', sv(k), Ts(k, i0), 10*log10(Ts(k, i0)), f(ip));
end
fprintf('quarter-wave thickness in a 2500 m/s layer: %.0f um at 5 MHz, %.1f um at 15 MHz\n', 2500/(4*5e6)*1e6, 2500/(4*15e6)*1e6);

%% 3. Layer attenuation (alpha in dB/mm/MHz, linear in f), loss at f0
c1 = 2500; al = [0 0.5 1 2];
for f0 = [5e6 15e6]
  d = c1/(4*f0);                                  % quarter-wave thickness (m)
  for a = al
    ad = a*(f0/1e6)*d*1e3/8.686;                  % Np through the layer
    M = lineM(Z1, pi/2, ad); v = M*[1; 1/ZL];   % [p_in; u_in] for p_L = 1, u_L = 1/ZL
    Zi = v(1)/v(2); absorb = (1/ZL)/real(v(1)*conj(v(2)));   % P_load / P_into_layer
    fprintf('f0 = %4.1f MHz, alpha = %.1f dB/mm/MHz: absorption %.2f dB, total one-way %.2f dB\n', f0/1e6, a, -10*log10(absorb), -10*log10(absorb*Tpow(Zi)));
  end
end

%% 4. Glue bondline between piezo and layer (Zb = 3 MRayl, cb = 2500 m/s)
f0s = linspace(1e6, 40e6, 400); bonds = [2 5 10]*1e-6; Tb = zeros(numel(bonds), numel(f0s));
for b = 1:numel(bonds)
  for k = 1:numel(f0s)
    Mb = lineM(3, 2*pi*f0s(k)*bonds(b)/2500, 0); Ml = lineM(Z1, pi/2, 0); M = Mb*Ml;
    Tb(b,k) = Tpow((M(1,1)*ZL + M(1,2))/(M(2,1)*ZL + M(2,2)));
  end
end
for fq = [5e6 20e6 30e6]
  [~, k] = min(abs(f0s - fq));
  fprintf('bondline at %2.0f MHz: 2 um %.2f dB, 5 um %.2f dB, 10 um %.2f dB\n', fq/1e6, 10*log10(Tb(:,k)));
end

%% Figures
figure('Position', [80 80 1000 330]);
subplot(1,3,1); plot(f, T1, 'LineWidth', 1.6); hold on; plot(f, T0 + 0*f, '--'); plot([band(1) band(end)], [0.5 0.5], 'k-o');
grid on; xlabel('f / f_0'); ylabel('power transmission'); title('One \lambda/4 layer, Z_1 = 7.16 MRayl'); legend('\lambda/4 layer', 'no layer', '-3 dB band', 'Location', 'south');
subplot(1,3,2); plot(f, Ts', 'LineWidth', 1.2); hold on; plot(f, T1, 'k', 'LineWidth', 1.6); grid on; xlim([0.5 1.5]); ylim([0.4 1.02]);
xlabel('f / f_0'); ylabel('power transmission'); title('Thickness error x0.9 ... x1.1'); legend('0.90', '0.95', '1.05', '1.10', 'nominal', 'Location', 'south');
subplot(1,3,3); plot(f0s/1e6, 10*log10(Tb'), 'LineWidth', 1.4); grid on; xlabel('design frequency f_0 (MHz)'); ylabel('transmission at f_0 (dB)');
title('Glue bondline (3 MRayl) under the layer'); legend('2 \mum', '5 \mum', '10 \mum', 'Location', 'southwest');
