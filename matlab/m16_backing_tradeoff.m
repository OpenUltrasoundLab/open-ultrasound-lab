%% m16_backing_tradeoff.m  Backing impedance vs bandwidth, insertion loss and ring-down (KLM/ABCD)
% Purpose : sweep the backing impedance ZB of a 5 MHz PZT-5H plate radiating into water, with no
%           matching layer, one layer (Zp^1/3 ZL^2/3) or two layers (Zp^4/7 ZL^3/7, Zp^1/7 ZL^6/7),
%           untuned, 50 ohm source and receiver. Report the -6 dB fractional bandwidth, the peak
%           insertion loss IL = 20 log10|2H| and the -40 dB envelope length, and check the mechanical
%           Q formula Q = (pi/2) Zp/(ZB + ZF) against the open-circuit KLM response.
% Change  : the layer set, material, R, the ZB grid.
% Expected: no matching: FBW 7 -> 19 -> 30 -> 58 -> 118 % and IL -5.6 -> -22.7 -> -31.2 -> -41.0 -> -48.1 dB
%           at ZB = 0, 5, 10, 20, 34.5 MRayl; two layers: 65-80 % at -17 to -34 dB.
clear; close all; rng(1);

e0   = 8.854187817e-12;
mm   = @(M,N) [M(:,1).*N(:,1)+M(:,2).*N(:,3), M(:,1).*N(:,2)+M(:,2).*N(:,4), ...
               M(:,3).*N(:,1)+M(:,4).*N(:,3), M(:,3).*N(:,2)+M(:,4).*N(:,4)];
ser  = @(Z) [ones(size(Z)), Z, zeros(size(Z)), ones(size(Z))];
shu  = @(Y) [ones(size(Y)), zeros(size(Y)), Y, ones(size(Y))];
xfm  = @(p) [p, zeros(size(p)), zeros(size(p)), 1./p];
tl   = @(Z,th) [cos(th), 1i*Z.*sin(th), 1i*sin(th)./Z, cos(th)];
zlin = @(Z,th,ZT) Z.*(ZT.*cos(th)+1i*Z.*sin(th))./(Z.*cos(th)+1i*ZT.*sin(th));
kof   = @(f,P) 2*pi*f ./ (P.v*sqrt(1+1i/P.Qm));
Z0of  = @(P) P.rho*P.v*sqrt(1+1i/P.Qm)*P.A;
C0of  = @(P) P.er*e0*(1-1i*P.tand)*P.A/P.d;
hof   = @(P) P.kt*sqrt(P.rho*P.v^2*P.er*e0)/(P.er*e0*(1-1i*P.tand));
X1of  = @(f,P) hof(P)^2*sin(kof(f,P)*P.d)./((2*pi*f).^2*Z0of(P));
phiof = @(f,P) 2*hof(P)*sin(kof(f,P)*P.d/2)./(2*pi*f*Z0of(P));
klm   = @(f,P,ZB) mm(mm(mm(mm(ser(1./(1i*2*pi*f*C0of(P))), ser(1i*X1of(f,P))), xfm(phiof(f,P))), ...
        shu(1./zlin(Z0of(P)*ones(size(f)), kof(f,P)*P.d/2, ZB*ones(size(f))))), tl(Z0of(P)*ones(size(f)), kof(f,P)*P.d/2));

fa = 5e6; P = struct('rho',7500,'v',4600,'kt',0.51,'er',1470,'Qm',65,'tand',0.02,'d',4600/(2*fa),'A',1.5e-6);
Zp = P.rho*P.v/1e6; ZLm = 1.5; ZL = ZLm*1e6*P.A; R = 50; fm = 0.9*fa;
sets = {'none', [], []; 'one layer', Zp^(1/3)*ZLm^(2/3), 2500; 'two layers', [Zp^(4/7)*ZLm^(3/7), Zp^(1/7)*ZLm^(6/7)], [2800 2400]};
Nfft = 1024; fs = 25.6*fa; df = fs/Nfft; fg = (0:Nfft/2)'*df; fg(1) = 1e-3*df; lim = round(2*fa/df);
wlp = ones(size(fg)); bnd = fg > 2*fa & fg < 2.5*fa; wlp(bnd) = 0.5*(1+cos(pi*(fg(bnd)-2*fa)/(0.5*fa))); wlp(fg >= 2.5*fa) = 0;
ZBs = 0:0.5:34.5;
res = zeros(numel(ZBs), 3, size(sets,1));                  % [FBW %, IL dB, -40 dB length us]

%% Sweep
for s = 1:size(sets,1)
    for q = 1:numel(ZBs)
        M = klm(fg, P, ZBs(q)*1e6*P.A);
        for L = 1:numel(sets{s,2})
            Zl = sets{s,2}(L); vl = sets{s,3}(L);
            M = mm(M, tl(Zl*1e6*P.A*ones(size(fg)), 2*pi*fg/vl*vl/(4*fm)));
        end
        Dn = M(:,1)*ZL + M(:,2) + R*(M(:,3)*ZL + M(:,4)); H = 2*R*ZL./Dn.^2;
        mag = abs(H(2:lim)); fb = fg(2:lim); [pk, ip] = max(mag); lev = pk/2;
        i1 = find(mag(1:ip) <= lev, 1, 'last'); i2 = ip - 1 + find(mag(ip:end) <= lev, 1, 'first');
        if isempty(i1), fl = fb(1); else fl = fb(i1) + (lev-mag(i1))/(mag(i1+1)-mag(i1))*(fb(i1+1)-fb(i1)); end
        fu = fb(i2-1) + (lev-mag(i2-1))/(mag(i2)-mag(i2-1))*(fb(i2)-fb(i2-1));
        Hs = zeros(Nfft,1); Hs(2:Nfft/2) = 2*H(2:Nfft/2).*wlp(2:Nfft/2); env = abs(ifft(Hs)); [epk, ie] = max(env);
        a = ie; while a > 1 && env(a) > 0.01*epk, a = a-1; end
        b = ie; while b < Nfft && env(b) > 0.01*epk, b = b+1; end
        res(q,:,s) = [100*(fu-fl)/((fu+fl)/2), 20*log10(2*pk), (b-a)/fs*1e6];
    end
end
fprintf('  ZB (MRayl) |  none: FBW  IL    t40  | one layer: FBW  IL    t40 | two layers: FBW  IL    t40\n');
for zb = [0 5 10 20 34.5]
    q = find(abs(ZBs-zb) < 1e-9);
    fprintf('   %5.1f     | %6.1f %6.1f %5.2f   | %8.1f %6.1f %5.2f     | %8.1f %6.1f %5.2f\n', zb, res(q,:,1), res(q,:,2), res(q,:,3));
end
fprintf('  (t40 = 8.00 us means the ring-down exceeds the %0.f us analysis window)\n', Nfft/fs*1e6);

%% Q check: open-circuit (R = 10 Mohm) transmit force of the bare, lossless plate
PL = P; PL.Qm = Inf; PL.tand = 0; fq = (3e6:2.5e3:7e6)';
ZBq = [0 1.5 3 5 8 12]; Qs = zeros(size(ZBq));
for q = 1:numel(ZBq)
    M = klm(fq, PL, ZBq(q)*1e6*P.A); D = M(:,1)*ZL + M(:,2) + 1e7*(M(:,3)*ZL + M(:,4)); F = abs(ZL./D);
    [pk, ip] = max(F); i1 = find(F(1:ip) < pk/sqrt(2), 1, 'last'); i2 = ip - 1 + find(F(ip:end) < pk/sqrt(2), 1, 'first');
    lv = pk/sqrt(2);   % interpolate the -3 dB crossings (grid samples alone bias Q low by ~1 %)
    fl = fq(i1) + (lv - F(i1))*(fq(i1+1) - fq(i1))/(F(i1+1) - F(i1));
    fu = fq(i2-1) + (lv - F(i2-1))*(fq(i2) - fq(i2-1))/(F(i2) - F(i2-1));
    Qs(q) = fq(ip)/(fu - fl);
end
Qth = pi/2*Zp./(ZBq + ZLm);
fprintf('\n  ZB    Q sim   Q = (pi/2) Zp/(ZB+ZF)   ratio\n');
fprintf('  %4.1f  %6.1f   %6.1f                %5.3f\n', [ZBq; Qs; Qth; Qs./Qth]);

%% Figures
figure('Name','backing trade-off');
subplot(2,1,1); plot(ZBs, res(:,1,1), 'k', ZBs, res(:,1,2), 'm', ZBs, res(:,1,3), 'b', 'LineWidth', 1.5); grid on;
xlabel('backing Z_B (MRayl)'); ylabel('-6 dB FBW (%)'); legend('no layer', 'one layer', 'two layers', 'Location', 'northwest'); title('Pulse-echo bandwidth');
subplot(2,1,2); plot(ZBs, res(:,2,1), 'k', ZBs, res(:,2,2), 'm', ZBs, res(:,2,3), 'b', 'LineWidth', 1.5); grid on;
xlabel('backing Z_B (MRayl)'); ylabel('IL (dB)'); title('Peak insertion loss, 50 \Omega, untuned');
figure('Name','Q formula');
zz = linspace(0, 12.5, 200);
plot(zz, pi/2*Zp./(zz+ZLm), 'b', 'LineWidth', 1.5); hold on; plot(ZBq, Qs, 'ro', 'MarkerFaceColor', 'r'); grid on;
xlabel('backing Z_B (MRayl)'); ylabel('mechanical Q'); legend('(\pi/2) Z_p/(Z_B + Z_F)', 'KLM, open circuit'); title('Q of a loaded half-wave plate (water front)');
