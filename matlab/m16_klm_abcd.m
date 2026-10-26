%% m16_klm_abcd.m  KLM transducer model with ABCD matrices: impedance, pulse-echo, impulse response
% Purpose : build the KLM model of a thickness-mode plate (Krimholtz, Leedom & Matthaei 1970) as a chain
%           of 2x2 ABCD matrices, add backing, two quarter-wave matching layers, a series tuning inductor
%           and a 50-ohm source/receiver, and compute the electrical input impedance, the pulse-echo
%           transfer function H = 2 R ZL / Delta^2 and the impulse response.
%           Four validation checks run first: free plate vs closed form, KLM vs Mason (independent nodal
%           solution), energy balance of a lossless chain, IEEE Std 176 kt recovery and the Q formula.
% Change  : material constants, fa, area A, backing ZB, layer set, tuning, R.
% Expected: validation errors ~1e-12 or smaller; reference design fc ~ 4.30 MHz, FBW ~ 67.5 %,
%           IL ~ -10.1 dB, envelope widths ~ 0.38 / 0.87 / 2.2 us at -6 / -20 / -40 dB.
% Conventions: acoustic "voltage" = force F (N), "current" = particle velocity u (m/s);
%           matrices are stored as Nf x 4 arrays [A B C D], one row per frequency.
clear; close all; rng(1);

%% ABCD building blocks (anonymous functions only: runs in MATLAB and Octave)
e0   = 8.854187817e-12;                                   % vacuum permittivity (F/m)
mm   = @(M,N) [M(:,1).*N(:,1)+M(:,2).*N(:,3), M(:,1).*N(:,2)+M(:,2).*N(:,4), ...
               M(:,3).*N(:,1)+M(:,4).*N(:,3), M(:,3).*N(:,2)+M(:,4).*N(:,4)];
ser  = @(Z) [ones(size(Z)), Z, zeros(size(Z)), ones(size(Z))];      % series impedance
shu  = @(Y) [ones(size(Y)), zeros(size(Y)), Y, ones(size(Y))];      % shunt admittance
xfm  = @(p) [p, zeros(size(p)), zeros(size(p)), 1./p];             % ideal transformer V = p F, I = u/p
tl   = @(Z,th) [cos(th), 1i*Z.*sin(th), 1i*sin(th)./Z, cos(th)];    % transmission line
zlin = @(Z,th,ZT) Z.*(ZT.*cos(th)+1i*Z.*sin(th))./(Z.*cos(th)+1i*ZT.*sin(th));   % line input impedance
zin  = @(M,ZL) (M(:,1).*ZL+M(:,2))./(M(:,3).*ZL+M(:,4));

%% Piezo plate: KLM elements from (rho, vD, kt, epsS_r, Qm, tan delta, thickness d, area A)
kof   = @(f,P) 2*pi*f ./ (P.v*sqrt(1+1i/P.Qm));                  % complex wavenumber (1/m)
Z0of  = @(P) P.rho*P.v*sqrt(1+1i/P.Qm)*P.A;                       % plate mechanical impedance (kg/s)
C0of  = @(P) P.er*e0*(1-1i*P.tand)*P.A/P.d;                       % clamped capacitance (F)
hof   = @(P) P.kt*sqrt(P.rho*P.v^2*P.er*e0)/(P.er*e0*(1-1i*P.tand));   % h33 = e33/epsS (V/m)
X1of  = @(f,P) hof(P)^2*sin(kof(f,P)*P.d)./((2*pi*f).^2*Z0of(P));     % KLM series reactance (ohm)
phiof = @(f,P) 2*hof(P)*sin(kof(f,P)*P.d/2)./(2*pi*f*Z0of(P));      % KLM transformer (V/N)
klm   = @(f,P,ZB) mm(mm(mm(mm(ser(1./(1i*2*pi*f*C0of(P))), ser(1i*X1of(f,P))), xfm(phiof(f,P))), ...
        shu(1./zlin(Z0of(P)*ones(size(f)), kof(f,P)*P.d/2, ZB*ones(size(f))))), tl(Z0of(P)*ones(size(f)), kof(f,P)*P.d/2));

fa = 5e6;                                                  % plate antiresonance vD/(2t) (Hz)
P  = struct('rho',7500,'v',4600,'kt',0.51,'er',1470,'Qm',65,'tand',0.02,'d',4600/(2*fa),'A',1.5e-6);  % PZT-5H, 0.3 x 5 mm
PL = P; PL.Qm = Inf; PL.tand = 0;                          % lossless copy for the analytic checks
Zp = P.rho*P.v/1e6;                                        % plate impedance (MRayl)
fprintf('PZT-5H plate: t = %.3f mm, Zp = %.1f MRayl, C0 = %.1f pF\n', P.d*1e3, Zp, real(C0of(P))*1e12);

%% Validation 1: free lossless plate vs closed form Z = (1/(j w C0)) (1 - kt^2 tan(x)/x)
f  = (0.2e6:13.7e3:20e6)';
Zk = zin(klm(f,PL,0), 0);
x  = 2*pi*f*PL.d/(2*PL.v);
Za = (1 - PL.kt^2*tan(x)./x) ./ (1i*2*pi*f*real(C0of(PL)));
fprintf('V1 free plate, KLM vs closed form   : max rel diff = %.1e\n', max(abs(Zk-Za)./abs(Za)));

%% Validation 2: KLM vs Mason (nodal analysis of the T-network + transformer N = h C0 and -C0)
Pm = P; Pm.Qm = 30;                                        % lossy, loaded plate
ZB = 5e6*P.A; ZF = 1.5e6*P.A;                              % backing and front loads (kg/s)
Zk = zin(klm(f,Pm,ZB), ZF);
k  = kof(f,Pm); Z0 = Z0of(Pm); C0 = C0of(Pm); N = hof(Pm)*C0;
ZT = 1i*Z0*tan(k*Pm.d/2); ZS = Z0./(1i*sin(k*Pm.d));      % Mason arms
Y  = 1./(ZT+ZB) + 1./(ZT+ZF); G = Y./(1+ZS.*Y);            % velocity into the T-network per transformer force
jwC = 1i*2*pi*f*C0;
Ip = N^2*G./(1 - N^2*G./jwC);                              % primary current per volt
Zm = 1./(jwC + Ip);
fprintf('V2 loaded lossy plate, KLM vs Mason : max rel diff = %.1e\n', max(abs(Zk-Zm)./abs(Zm)));

%% Reference stack: backing, two Desilets-type layers at fm = 0.9 fa, series L, R = 50 ohm, water load
ZBm = 5; ZLm = 1.5; R = 50; fm = 0.9*fa;
Zl  = [Zp^(4/7)*ZLm^(3/7), Zp^(1/7)*ZLm^(6/7)];           % layer impedances (MRayl): 9.00, 2.35
vl  = [2800, 2400]; tlay = vl/(4*fm);                      % layer speeds (m/s), quarter-wave thicknesses (m)
Lser = 1/((2*pi*fm)^2*real(C0of(P)));                     % series tuning inductor (H)
stack = @(f,Pp,ZBk,Lt) mm(ser(1i*2*pi*f*Lt), mm(mm(klm(f,Pp,ZBk*1e6*Pp.A), tl(Zl(1)*1e6*Pp.A*ones(size(f)), 2*pi*f/vl(1)*tlay(1))), ...
                     tl(Zl(2)*1e6*Pp.A*ones(size(f)), 2*pi*f/vl(2)*tlay(2))));
fprintf('Layers: Z = %.2f, %.2f MRayl, t = %.3f, %.3f mm; series L = %.1f uH\n', Zl, tlay*1e3, Lser*1e6);

%% Validation 3: energy balance of the lossless chain: Re{V I*} = Re{ZL} |u|^2 with u = 1
fe = fa*(0.05 + 1.8*(0:159)'/160) + 917;
M  = stack(fe, PL, 0, Lser); ZL = ZLm*1e6*P.A;
V  = M(:,1)*ZL + M(:,2); I = M(:,3)*ZL + M(:,4);
fprintf('V3 energy balance, lossless chain   : max rel diff = %.1e\n', max(abs(real(V.*conj(I)) - real(ZL))/real(ZL)));

%% Validation 4: IEEE Std 176 on the simulated free plate, and the open-circuit Q formula
imZ = @(ff) imag(zin(klm(ff,PL,0),0)); imY = @(ff) imag(1./zin(klm(ff,PL,0),0));
a = 0.3*fa; b = 0.999*fa; for it = 1:60, c = (a+b)/2; if imZ(a)*imZ(c) <= 0, b = c; else a = c; end, end, fr = (a+b)/2;
a = 1.001*fr; b = 1.5*fa;  for it = 1:60, c = (a+b)/2; if imY(a)*imY(c) <= 0, b = c; else a = c; end, end, fan = (a+b)/2;
kt176 = sqrt(pi/2*fr/fan*tan(pi/2*(fan-fr)/fan));
fprintf('V4 fr = %.4f MHz, fa = %.4f MHz, IEEE 176 kt = %.4f (input %.2f)\n', fr/1e6, fan/1e6, kt176, PL.kt);
fq = (3e6:2.5e3:7e6)';
for ZBq = [0 5]
    Mq = klm(fq, PL, ZBq*1e6*P.A); Dq = Mq(:,1)*ZL + Mq(:,2) + 1e7*(Mq(:,3)*ZL + Mq(:,4));
    Ft = abs(ZL./Dq); [pk, ip] = max(Ft); lev = pk/sqrt(2);
    i1 = find(Ft(1:ip) < lev, 1, 'last'); i2 = ip - 1 + find(Ft(ip:end) < lev, 1, 'first');
    fl = fq(i1) + (lev - Ft(i1))*(fq(i1+1) - fq(i1))/(Ft(i1+1) - Ft(i1));   % interpolated -3 dB crossings
    fu = fq(i2-1) + (lev - Ft(i2-1))*(fq(i2) - fq(i2-1))/(Ft(i2) - Ft(i2-1));
    Qs = fq(ip)/(fu - fl);
    fprintf('   open-circuit Q, ZB = %g MRayl: simulated %.2f, (pi/2) Zp/(ZB+ZF) = %.2f\n', ZBq, Qs, pi/2*Zp/(ZBq+ZLm));
end

%% Reference design: spectrum on a grid scaled to fa, -6 dB band, impulse response
Nfft = 1024; fs = 25.6*fa; df = fs/Nfft; fg = (0:Nfft/2)'*df; fg(1) = 1e-3*df;
M  = stack(fg, P, ZBm, Lser);
Dn = M(:,1)*ZL + M(:,2) + R*(M(:,3)*ZL + M(:,4));
H  = 2*R*ZL./Dn.^2;                                        % pulse-echo V_out/V_source, unit planar reflector
Zi = zin(M, ZL);
lim = round(2*fa/df); mag = abs(H(2:lim)); fb = fg(2:lim);
[pk, ip] = max(mag); lev = pk/2;
i1 = find(mag(1:ip) <= lev, 1, 'last'); i2 = ip - 1 + find(mag(ip:end) <= lev, 1, 'first');
fl = fb(i1) + (lev-mag(i1))/(mag(i1+1)-mag(i1))*(fb(i1+1)-fb(i1));
fu = fb(i2-1) + (lev-mag(i2-1))/(mag(i2)-mag(i2-1))*(fb(i2)-fb(i2-1));
fc = (fl+fu)/2; fbw = (fu-fl)/fc; IL = 20*log10(2*pk);
wlp = ones(size(fg)); band2 = fg > 2*fa & fg < 2.5*fa; wlp(band2) = 0.5*(1+cos(pi*(fg(band2)-2*fa)/(0.5*fa))); wlp(fg >= 2.5*fa) = 0;
Hs = zeros(Nfft,1); Hs(2:Nfft/2) = 2*H(2:Nfft/2).*wlp(2:Nfft/2);    % analytic spectrum, low-passed 2-2.5 fa
asig = ifft(Hs); env = abs(asig); ir = real(asig); [epk, ie] = max(env);
dur = zeros(1,3); lv = [0.5 0.1 0.01];
for q = 1:3
    a = ie; while a > 1 && env(a) > epk*lv(q), a = a-1; end
    b = ie; while b < Nfft && env(b) > epk*lv(q), b = b+1; end
    dur(q) = (b-a)/fs;
end
fprintf('Reference design: fc = %.2f MHz, -6 dB FBW = %.1f %%, IL = %.1f dB at %.2f MHz\n', fc/1e6, 100*fbw, IL, fb(ip)/1e6);
fprintf('   envelope widths -6/-20/-40 dB = %.2f / %.2f / %.2f us\n', dur*1e6);
ifc = round(fc/df)+1;
fprintf('   Zin(fc) = %.0f %+.0fj ohm (untuned element alone at fm: ', real(Zi(ifc)), imag(Zi(ifc)));
Mu = stack(fm, P, ZBm, 0); zu = zin(Mu, ZL); fprintf('%.0f %+.0fj ohm)\n', real(zu), imag(zu));

%% Figures
figure('Name','KLM: impedance');
sel = fg <= 2.2*fa & fg > 0.05*fa;
subplot(2,1,1); semilogy(fg(sel)/1e6, abs(Zi(sel)), 'b', 'LineWidth', 1.5); grid on;
xlabel('frequency (MHz)'); ylabel('|Z_{in}| (\Omega)'); title('Electrical input impedance, reference stack (series L)');
subplot(2,1,2); plot(fg(sel)/1e6, angle(Zi(sel))*180/pi, 'r', 'LineWidth', 1.2); grid on;
xlabel('frequency (MHz)'); ylabel('phase (deg)');
figure('Name','KLM: pulse-echo');
sel = fg <= 3.5*fa & fg > 0;
subplot(2,1,1); plot(fg(sel)/1e6, 20*log10(2*abs(H(sel))), 'b', 'LineWidth', 1.5); grid on; hold on;
plot([fl fl]/1e6, [-80 0], 'k--', [fu fu]/1e6, [-80 0], 'k--'); ylim([IL-60, IL+5]);
xlabel('frequency (MHz)'); ylabel('20 log_{10}|2H| (dB)'); title(sprintf('Pulse-echo: f_c = %.2f MHz, FBW = %.1f %%, IL = %.1f dB', fc/1e6, 100*fbw, IL));
idx = mod((ie-60:ie+300)-1, Nfft) + 1; tt = ((-60:300)/fs)*1e6;
subplot(2,1,2); plot(tt, ir(idx)/epk, 'b', tt, env(idx)/epk, 'r', 'LineWidth', 1.2); grid on;
xlabel('time re envelope peak (\mus)'); ylabel('normalized'); title('Pulse-echo impulse response (low-passed at 2-2.5 f_a)');
