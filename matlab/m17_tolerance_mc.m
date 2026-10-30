%% m17_tolerance_mc.m  --  Monte Carlo tolerance study of a two-layer matched stack
% Purpose : build N virtual transducers with random layer thickness, layer
%           impedance, bond-line thickness, plate thickness and kerf (element
%           area) and read the spread of centre frequency, -6 dB fractional
%           bandwidth and two-way insertion loss from the KLM + ABCD chain
%           of m17_matching_layers.m (Section 17.10, Demo 17.8).
% Model   : PZT-5H-like plate (lambda/2 at 5 MHz), reference port 20 mm^2
%           scaled by element width (pitch 0.30 mm, kerf 25 um nominal),
%           50 ohm, backing 5 MRayl, load 1.5 MRayl, two Desilets layers
%           (8.95 / 2.34 MRayl) at 2500 m/s, a 3 MRayl / 2500 m/s bond line
%           under each layer (truncated at zero), lossless layers.
% Change  : the sigma values below, Ntr, the window limits.
% Expect  : design 4.67 MHz / 66 % / -5.8 dB; built mean about 4.59 MHz
%           (bond lines shift it down by about 1.7 %), sigma about 0.06 MHz;
%           tornado: bond lines, then layer and plate thickness dominate fc.
%           Runtime ~ 13 s (Octave).
clear; close all; rng(1);

%% Parameters
rho = 7500; v = 4560; kt = 0.50; epsS = 1470*8.854e-12;
f0 = 5e6; t0 = v/(2*f0); A0 = 20e-6; ZB = 5; ZL = 1.5; Rs = 50; cLay = 2500;
Zp = rho*v/1e6; Zd = [Zp^(4/7)*ZL^(3/7), Zp^(1/7)*ZL^(6/7)];   % Desilets, MRayl
f = 1.5e6:0.05e6:8.5e6;                                      % scan grid (Hz)
Ntr = 240;                                                   % virtual builds
sd = 0.03; sz = 0.05; bm = 3e-6; bs = 1.5e-6; st = 0.01; sk = 3e-6;  % tolerances

lineM = @(Z, d, c, fr) [cos(2*pi*fr*d/c), 1i*Z*sin(2*pi*fr*d/c); 1i*sin(2*pi*fr*d/c)/Z, cos(2*pi*fr*d/c)];
zin = @(Mx, Zt) (Mx(1,1)*Zt + Mx(1,2)) / (Mx(2,1)*Zt + Mx(2,2));

% each row of allP = [dt kerf b1 b2 z1 z2 d1 d2] -> out = [fc (MHz), FBW (%), IL (dB)]
nom = [0 25e-6 0 0 0 0 0 0];
P = zeros(Ntr, 8);
for i = 1:Ntr
  P(i,:) = [st*randn, 25e-6 + sk*randn, max(0, bm + bs*randn), max(0, bm + bs*randn), sz*randn, sz*randn, sd*randn, sd*randn];
end
base = [0 25e-6 bm bm 0 0 0 0];
tor = repmat(base, 5, 1);
tor(1,7:8) = sd; tor(2,5:6) = sz; tor(3,3:4) = bm + bs; tor(4,1) = st; tor(5,2) = 25e-6 + sk;
allP = [nom; base; P; tor];
out = zeros(size(allP,1), 3);
for q = 1:size(allP,1)
  p = allP(q,:); t = t0*(1 + p(1)); A = A0*(0.300e-3 - p(2))/0.275e-3;
  Ls = zeros(0,3);                                           % rows: [Z (MRayl), d (m), c]
  if p(3) > 0, Ls(end+1,:) = [3, p(3), 2500]; end
  Ls(end+1,:) = [Zd(1)*(1 + p(5)), (1 + p(7))*cLay/(4*f0), cLay];
  if p(4) > 0, Ls(end+1,:) = [3, p(4), 2500]; end
  Ls(end+1,:) = [Zd(2)*(1 + p(6)), (1 + p(8))*cLay/(4*f0), cLay];
  C0 = epsS*A/t; Z0 = rho*v*A; h2 = kt^2*rho*v^2/epsS; m = zeros(size(f));
  for k = 1:numel(f)
    w = 2*pi*f(k); th = w*t/v; X1 = h2*sin(th)/(w^2*Z0); phi = (w*Z0/(2*sqrt(h2)))/sin(th/2);
    half = lineM(Z0, t/2, v, f(k)); Zb = zin(half, ZB*1e6*A); Mf = half;
    for L = 1:size(Ls,1), Mf = Mf*lineM(Ls(L,1)*1e6*A, Ls(L,2), Ls(L,3), f(k)); end
    T = [1, 1i*(X1 - 1/(w*C0)); 0, 1] * [1/phi, 0; 0, phi] * [1, 0; 1/Zb, 1] * Mf;
    ZLa = ZL*1e6*A; u = 1/((T(1,1)*ZLa + T(1,2)) + Rs*(T(2,1)*ZLa + T(2,2)));
    m(k) = abs(2*Rs/ZLa*(ZLa*u)^2);
  end
  [pk, ip] = max(m); lo = ip; while lo > 1 && m(lo) > pk/2, lo = lo - 1; end
  hi = ip; while hi < numel(m) && m(hi) > pk/2, hi = hi + 1; end
  fl = interp1(m(lo:lo+1), f(lo:lo+1), pk/2); fh = interp1(m(hi-1:hi), f(hi-1:hi), pk/2);
  out(q,:) = [(fl + fh)/2/1e6, 100*(fh - fl)/((fl + fh)/2), 20*log10(2*pk)];
end
D = out(1,:); M = out(2,:); B = out(3:2+Ntr,:); TT = out(3+Ntr:end,:);

%% Results
fprintf('design (no bond lines): fc %.2f MHz, FBW %.0f %%, IL %.1f dB\n', D);
fprintf('mean build (3 um bonds): fc %.2f MHz, FBW %.0f %%, IL %.1f dB\n', M);
fprintf('Monte Carlo (%d builds): fc %.2f +- %.3f MHz, FBW %.1f +- %.1f %%, IL %.2f +- %.2f dB\n', ...
        Ntr, mean(B(:,1)), std(B(:,1)), mean(B(:,2)), std(B(:,2)), mean(B(:,3)), std(B(:,3)));
win = B(:,1) >= 0.95*D(1) & B(:,1) <= 1.05*D(1) & B(:,2) >= D(2) - 5 & abs(B(:,3) - D(3)) <= 1.5;
fprintf('inside the three windows: %.0f %%\n', 100*mean(win));
names = {'layer thickness', 'layer Z', 'bond lines', 'plate thickness', 'kerf'};
for k = 1:5, fprintf('+1 sigma %-16s: dfc %+6.1f kHz, dIL %+5.2f dB\n', names{k}, 1000*(TT(k,1) - M(1)), TT(k,3) - M(3)); end

%% Figure
figure('Position', [100 100 900 320]);
subplot(1,3,1); hist(B(:,1), 20); hold on; plot([D(1) D(1)], ylim, 'k--'); xlabel('centre frequency (MHz)'); ylabel('builds'); title('f_c');
subplot(1,3,2); hist(B(:,2), 20); hold on; plot([D(2) D(2)], ylim, 'k--'); xlabel('-6 dB FBW (%)'); title('bandwidth');
subplot(1,3,3); barh(1000*(TT(:,1) - M(1))); set(gca, 'YTickLabel', names); xlabel('f_c shift per +1\sigma (kHz)'); title('tornado');
