%% m17_matching_layers.m  --  KLM + ABCD model of a backed plate with 0-3 matching layers
% Purpose : compare matching-layer design rules (quarter-wave sqrt, Desilets
%           1978, binomial) with one self-contained KLM/ABCD model and print
%           centre frequency, -6 dB fractional bandwidth, -20 dB pulse length
%           and two-way insertion loss (perfect reflector at the face).
% Model   : PZT-5H-like plate (rho 7500 kg/m^3, v 4560 m/s, kt 0.50,
%           eps_S/eps0 1470), lambda/2 at 5 MHz, area 20 mm^2, 50 ohm
%           source and receiver, backing 5 MRayl, load 1.5 MRayl, layers
%           lambda/4 at 5 MHz in a c = 2500 m/s material, lossless.
% Change  : ZB (backing), Rs, the design list, layer loss alpha_dB.
% Expect  : none 27 % / 1.25 us / -12.8 dB; sqrt 56 % / 1.03 us / -5.6 dB;
%           Desilets-1 44 % / 0.78 us / -5.4 dB; Desilets-2 66 % / 0.81 us;
%           binomial-2 79 % / 1.00 us; shortest single-layer pulse at
%           4.25-4.5 MRayl (the Desilets value). Runtime ~ 8 s (Octave).
% Note    : fc, FBW and loss are read in the fundamental band (peak searched
%           below 9 MHz), so the 3rd-harmonic resonance near 15 MHz is never
%           mistaken for the main band (it matters for light backings).
clear; close all; rng(1);

%% Parameters
rho = 7500; v = 4560; kt = 0.50; epsS = 1470*8.854e-12;    % piezo
f0 = 5e6; t = v/(2*f0); A = 20e-6;                         % plate thickness (m), area (m^2)
Zp = rho*v/1e6;                                            % 34.2 MRayl
ZB = 5; ZL = 1.5; Rs = 50;                                 % MRayl, MRayl, ohm
cLay = 2500; alpha_dB = 0;                                 % layer speed (m/s), loss (dB/mm/MHz)
N = 2048; fs = 100e6; df = fs/N; f = (1:N/2-1)*df;         % frequency grid (Hz)

% one ABCD for a (lossy) line: Z (Rayl*m^2), length d (m), speed c, loss a (Np/m)
lineM = @(Z, d, c, fr, a) [cosh((a + 1i*2*pi*fr/c)*d), Z*sinh((a + 1i*2*pi*fr/c)*d); ...
                           sinh((a + 1i*2*pi*fr/c)*d)/Z, cosh((a + 1i*2*pi*fr/c)*d)];
zin = @(Mx, Zt) (Mx(1,1)*Zt + Mx(1,2)) / (Mx(2,1)*Zt + Mx(2,2));
% -20 dB envelope duration with linear interpolation of both crossings (us)
dur20 = @(e, i1, i2, th) ((i2 + (e(i2) - th)/(e(i2) - e(i2+1))) - (i1 - 1 + (th - e(i1-1))/(e(i1) - e(i1-1))))/fs*1e6;

%% Designs (impedances in MRayl)
R = @(e) Zp^e * ZL^(1-e);
names = {'none', 'sqrt (QW)', 'Desilets-1', 'Desilets-2', 'binomial-2', 'binomial-3', 'grid best-2'};
Zs = {[], R(1/2), R(1/3), [R(4/7) R(1/7)], [R(3/4) R(1/4)], [R(7/8) R(1/2) R(1/8)], [8.0 2.4]};

C0 = epsS*A/t; Z0 = rho*v*A; h2 = kt^2*rho*v^2/epsS;   % clamped capacitance, line impedance, h^2
Hpe = zeros(numel(Zs), N/2-1); eta = Hpe;
for q = 1:numel(Zs)
  for k = 1:numel(f)
    w = 2*pi*f(k); th = w*t/v;
    X1 = h2*sin(th)/(w^2*Z0); phi = (w*Z0/(2*sqrt(h2)))/sin(th/2);
    half = lineM(Z0, t/2, v, f(k), 0);
    Zb = zin(half, ZB*1e6*A);                          % back branch seen from the centre node
    Mf = half;
    for L = 1:numel(Zs{q})
      Mf = Mf*lineM(Zs{q}(L)*1e6*A, cLay/(4*f0), cLay, f(k), alpha_dB*f(k)/1e6*1000/8.686);
    end
    T = [1, 1i*(X1 - 1/(w*C0)); 0, 1] * [1/phi, 0; 0, phi] * [1, 0; 1/Zb, 1] * Mf;
    ZLa = ZL*1e6*A;
    u = 1/((T(1,1)*ZLa + T(1,2)) + Rs*(T(2,1)*ZLa + T(2,2)));   % load velocity per source volt
    F = ZLa*u;
    eta(q,k) = 4*Rs*abs(u)^2*ZLa;                     % one-way transducer power gain
    Hpe(q,k) = 2*Rs/ZLa*F^2;                           % pulse-echo Vout/Vs (reciprocity, det = 1)
  end
end

%% Metrics
tvec = (0:N-1)/fs*1e6; tau0 = 2e-6;                    % time axis (us), display delay (s)
% receive low-pass for the time-domain view: 1 below 9 MHz, cos^2 taper to 0 at 12 MHz
% (removes the 3rd-harmonic thickness resonance near 14 MHz; the fundamental band is untouched)
Glp = (f < 9e6) + (f >= 9e6 & f < 12e6).*cos(pi/2*(f - 9e6)/3e6).^2;
fprintf('%-12s %-20s %7s %7s %9s %8s\n', 'design', 'Z (MRayl)', 'fc MHz', 'FBW %', 'L-20 us', 'IL2 dB');
rf = zeros(numel(Zs), N);
for q = 1:numel(Zs)
  m = abs(Hpe(q,:)); [pk, ip] = max(m.*(f < 9e6));             % fundamental band only
  lo = ip; while lo > 1 && m(lo) > pk/2, lo = lo - 1; end
  hi = ip; while hi < numel(m) && m(hi) > pk/2, hi = hi + 1; end
  fl = interp1(m(lo:lo+1), f(lo:lo+1), pk/2); fh = interp1(m(hi-1:hi), f(hi-1:hi), pk/2);
  fc = (fl + fh)/2; fbw = (fh - fl)/fc;
  S = zeros(1, N); S(2:N/2) = 2*Hpe(q,:).*Glp.*exp(-1i*2*pi*f*tau0);   % analytic-signal spectrum
  a = ifft(S); env = abs(a); rf(q,:) = real(a);
  idx = find(env > 0.1*max(env)); L20 = dur20(env, idx(1), idx(end), 0.1*max(env));
  zstr = sprintf('%.2f ', Zs{q}); if isempty(Zs{q}), zstr = '-'; end
  fprintf('%-12s %-20s %7.2f %7.0f %9.2f %8.1f\n', names{q}, zstr, fc/1e6, 100*fbw, L20, 20*log10(2*pk));
end

%% Single layer: which Z gives the shortest pulse?
Zgrid = 2:0.25:12; L20g = zeros(size(Zgrid)); fbwg = L20g;
for g = 1:numel(Zgrid)
  H = zeros(1, numel(f));
  for k = 1:numel(f)
    w = 2*pi*f(k); th = w*t/v; X1 = h2*sin(th)/(w^2*Z0); phi = (w*Z0/(2*sqrt(h2)))/sin(th/2);
    half = lineM(Z0, t/2, v, f(k), 0); Zb = zin(half, ZB*1e6*A);
    Mf = half*lineM(Zgrid(g)*1e6*A, cLay/(4*f0), cLay, f(k), 0);
    T = [1, 1i*(X1 - 1/(w*C0)); 0, 1] * [1/phi, 0; 0, phi] * [1, 0; 1/Zb, 1] * Mf;
    ZLa = ZL*1e6*A; u = 1/((T(1,1)*ZLa + T(1,2)) + Rs*(T(2,1)*ZLa + T(2,2)));
    H(k) = 2*Rs/ZLa*(ZLa*u)^2;
  end
  S = zeros(1, N); S(2:N/2) = 2*H.*Glp.*exp(-1i*2*pi*f*tau0); env = abs(ifft(S));
  idx = find(env > 0.1*max(env)); L20g(g) = dur20(env, idx(1), idx(end), 0.1*max(env));
  m = abs(H); [pk, ip] = max(m.*(f < 9e6)); lo = ip; while lo > 1 && m(lo) > pk/2, lo = lo - 1; end
  hi = ip; while hi < numel(m) && m(hi) > pk/2, hi = hi + 1; end
  fbwg(g) = (f(hi) - f(lo))/((f(hi) + f(lo))/2);
end
[Lmin, ib] = min(L20g);
fprintf('\nSingle layer: shortest -20 dB pulse %.2f us at Z1 = %.2f MRayl (Desilets %.2f, sqrt %.2f)\n', Lmin, Zgrid(ib), R(1/3), R(1/2));
fprintf('  pulse length at Z1 = 4.25: %.2f us, at 7.25: %.2f us (the minimum is flat from 4 to 5 MRayl)\n', L20g(Zgrid == 4.25), L20g(Zgrid == 7.25));
[~, iw] = max(fbwg);
fprintf('Single layer: widest -6 dB band %.0f %% at Z1 = %.2f MRayl (ripple counts as band)\n', 100*fbwg(iw), Zgrid(iw));

%% Figures
figure('Position', [80 80 980 640]);
subplot(2,2,1); plot(f/1e6, 20*log10(eta([1 2 3 4],:))', 'LineWidth', 1.3); grid on; xlim([0 10]); ylim([-40 0]);
xlabel('frequency (MHz)'); ylabel('20 log_{10} \eta(f)  (dB, two-way)'); title('Pulse-echo loss, backing 5 MRayl');
legend(names([1 2 3 4]), 'Location', 'south');
subplot(2,2,2); plot(f/1e6, 20*log10(eta([4 5 6 7],:))', 'LineWidth', 1.3); grid on; xlim([0 10]); ylim([-40 0]);
xlabel('frequency (MHz)'); ylabel('dB'); title('Two and three layers'); legend(names([4 5 6 7]), 'Location', 'south');
subplot(2,2,3); sel = [1 3 4 5];
for k = 1:4, plot(tvec - tau0*1e6, rf(sel(k),:)/max(abs(rf(sel(k),:))) - 2.2*(k-1), 'LineWidth', 1.1); hold on; end
xlim([-0.3 3]); grid on; xlabel('time (\mus)'); ylabel('normalised, offset'); title('Pulse-echo impulse responses');
legend(names(sel), 'Location', 'southeast');
subplot(2,2,4); [ax, h1, h2] = plotyy(Zgrid, L20g, Zgrid, 100*fbwg); grid on;
xlabel('single-layer Z_1 (MRayl)'); ylabel(ax(1), '-20 dB pulse length (\mus)'); ylabel(ax(2), '-6 dB FBW (%)');
title('One layer: shortest pulse vs widest band');
