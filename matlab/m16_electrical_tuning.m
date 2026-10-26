%% m16_electrical_tuning.m  Tuning inductors, a transformer and a coaxial cable as ABCD factors
% Purpose : take the reference stack (PZT-5H plate, fa = 5 MHz, 1.5 mm^2, ZB = 5 MRayl, two layers at
%           4.5 MHz, water load) and compare electrical networks between it and the system:
%           none, series L, shunt L, series L + ideal transformer, and a 50-ohm coaxial cable
%           (100 pF/m) into a 50-ohm or a 1-kohm receiver. Prints peak two-way gain 20 log10|2 Vout/Vs|
%           (equal to the insertion loss when the source and receiver impedances are equal).
% Change  : transformer ratio n, cable length, receiver impedance, which capacitance L resonates.
% Expected: untuned -21.6 dB, series L -10.1 dB, shunt L within 0.03 dB of untuned (50 ohm source),
%           series L + n = 0.4 about -2.6 dB; 50-ohm system: cable length has no effect on |H|;
%           1-kohm receiver: 2 m of cable costs about 9.6 dB.
clear; close all; rng(1);

e0   = 8.854187817e-12;
mm   = @(M,N) [M(:,1).*N(:,1)+M(:,2).*N(:,3), M(:,1).*N(:,2)+M(:,2).*N(:,4), ...
               M(:,3).*N(:,1)+M(:,4).*N(:,3), M(:,3).*N(:,2)+M(:,4).*N(:,4)];
ser  = @(Z) [ones(size(Z)), Z, zeros(size(Z)), ones(size(Z))];
shu  = @(Y) [ones(size(Y)), zeros(size(Y)), Y, ones(size(Y))];
xfm  = @(p) [p, zeros(size(p)), zeros(size(p)), 1./p];
tl   = @(Z,th) [cos(th), 1i*Z.*sin(th), 1i*sin(th)./Z, cos(th)];
zlin = @(Z,th,ZT) Z.*(ZT.*cos(th)+1i*Z.*sin(th))./(Z.*cos(th)+1i*ZT.*sin(th));
zin  = @(M,ZL) (M(:,1).*ZL+M(:,2))./(M(:,3).*ZL+M(:,4));
kof   = @(f,P) 2*pi*f ./ (P.v*sqrt(1+1i/P.Qm));
Z0of  = @(P) P.rho*P.v*sqrt(1+1i/P.Qm)*P.A;
C0of  = @(P) P.er*e0*(1-1i*P.tand)*P.A/P.d;
hof   = @(P) P.kt*sqrt(P.rho*P.v^2*P.er*e0)/(P.er*e0*(1-1i*P.tand));
X1of  = @(f,P) hof(P)^2*sin(kof(f,P)*P.d)./((2*pi*f).^2*Z0of(P));
phiof = @(f,P) 2*hof(P)*sin(kof(f,P)*P.d/2)./(2*pi*f*Z0of(P));
klm   = @(f,P,ZB) mm(mm(mm(mm(ser(1./(1i*2*pi*f*C0of(P))), ser(1i*X1of(f,P))), xfm(phiof(f,P))), ...
        shu(1./zlin(Z0of(P)*ones(size(f)), kof(f,P)*P.d/2, ZB*ones(size(f))))), tl(Z0of(P)*ones(size(f)), kof(f,P)*P.d/2));

fa = 5e6; P = struct('rho',7500,'v',4600,'kt',0.51,'er',1470,'Qm',65,'tand',0.02,'d',4600/(2*fa),'A',1.5e-6);
Zp = P.rho*P.v/1e6; ZLm = 1.5; ZL = ZLm*1e6*P.A; Rs = 50; fm = 0.9*fa;
Zl = [Zp^(4/7)*ZLm^(3/7), Zp^(1/7)*ZLm^(6/7)]; vl = [2800 2400];
Nfft = 1024; fs = 25.6*fa; df = fs/Nfft; f = (0:Nfft/2)'*df; f(1) = 1e-3*df; lim = round(2*fa/df);
Mt = mm(mm(klm(f, P, 5e6*P.A), tl(Zl(1)*1e6*P.A*ones(size(f)), 2*pi*f/vl(1)*vl(1)/(4*fm))), tl(Zl(2)*1e6*P.A*ones(size(f)), 2*pi*f/vl(2)*vl(2)/(4*fm)));
C0 = real(C0of(P)); Lopt = 1/((2*pi*fm)^2*C0); w = 2*pi*f; I4 = [ones(size(f)) zeros(size(f)) zeros(size(f)) ones(size(f))];
cable = @(len) tl(50*ones(size(f)), w*len/2e8);            % 50-ohm coax, 2e8 m/s -> 100 pF/m
fprintf('C0 = %.1f pF, Lopt(4.5 MHz) = %.1f uH, cable C = 100 pF/m\n', C0*1e12, Lopt*1e6);

% gain(E, Rr): two-way voltage gain spectrum with electrical network E in front of the transducer
Dfun = @(M,R) M(:,1)*ZL + M(:,2) + R*(M(:,3)*ZL + M(:,4));
gain = @(E,Rr) (ZL./Dfun(mm(E,Mt),Rs)) .* (2*Rr./Dfun(mm(E,Mt),Rr));
peak = @(H) 20*log10(2*max(abs(H(2:lim))));

%% Cases
Hn  = gain(I4, Rs);
Hse = gain(ser(1i*w*Lopt), Rs);
Hsh = gain(shu(1./(1i*w*Lopt)), Rs);
Hx  = gain(mm(xfm(0.4*ones(size(f))), ser(1i*w*Lopt)), Rs);
fprintf('\n 50-ohm system                          peak gain (dB)\n');
fprintf('  untuned                               %7.2f\n', peak(Hn));
fprintf('  series L                              %7.2f\n', peak(Hse));
fprintf('  shunt L                               %7.2f  (difference from untuned %.4f dB)\n', peak(Hsh), peak(Hsh)-peak(Hn));
for n = [0.5 0.4 0.35]
    fprintf('  series L + transformer n = %.2f      %7.2f\n', n, peak(gain(mm(xfm(n*ones(size(f))), ser(1i*w*Lopt)), Rs)));
end
fprintf('\n cable length (m)   50-ohm rx, series L   1-kohm rx, untuned   1-kohm rx, shunt L(C0+Ccab)\n');
for len = [0 1 2 3]
    Lsh = 1/((2*pi*fm)^2*(C0 + len*100e-12));
    fprintf('   %3.1f              %7.2f               %7.2f              %7.2f\n', len, peak(gain(mm(cable(len), ser(1i*w*Lopt)), Rs)), ...
        peak(gain(cable(len), 1000)), peak(gain(mm(cable(len), shu(1./(1i*w*Lsh))), 1000)));
end
fprintf(' lumped open-circuit estimate for 2 m: 20 log10(C0/(C0+200 pF)) = %.1f dB\n', 20*log10(C0/(C0+200e-12)));
Zu = zin(Mt, ZL); Zs = zin(mm(ser(1i*w*Lopt), Mt), ZL); Zx = zin(mm(mm(xfm(0.4*ones(size(f))), ser(1i*w*Lopt)), Mt), ZL);
i45 = round(fm/df) + 1;
fprintf('\n Z seen at 4.5 MHz: untuned %.0f%+.0fj, series L %.0f%+.0fj, + n = 0.4 %.0f%+.0fj ohm\n', real(Zu(i45)), imag(Zu(i45)), real(Zs(i45)), imag(Zs(i45)), real(Zx(i45)), imag(Zx(i45)));

%% Figures
sel = f > 0 & f <= 10e6;
figure('Name','tuning: gain');
plot(f(sel)/1e6, 20*log10(2*abs(Hn(sel))), 'k', f(sel)/1e6, 20*log10(2*abs(Hse(sel))), 'b', f(sel)/1e6, 20*log10(2*abs(Hx(sel))), 'g', f(sel)/1e6, 20*log10(2*abs(Hsh(sel))), 'r--', 'LineWidth', 1.4);
grid on; ylim([-60 0]); xlabel('frequency (MHz)'); ylabel('20 log_{10}|2V_{out}/V_s| (dB)');
legend('untuned', 'series L', 'series L + n = 0.4', 'shunt L', 'Location', 'south'); title('Electrical tuning of the reference element, 50 \Omega');
figure('Name','tuning: impedance seen');
subplot(2,1,1); plot(f(sel)/1e6, real(Zu(sel)), 'k', f(sel)/1e6, real(Zs(sel)), 'b', f(sel)/1e6, real(Zx(sel)), 'g', 'LineWidth', 1.3); grid on; ylim([0 1500]);
ylabel('Re Z (\Omega)'); title('Impedance seen by the system'); legend('untuned', 'series L', '+ n = 0.4');
subplot(2,1,2); plot(f(sel)/1e6, imag(Zu(sel)), 'k', f(sel)/1e6, imag(Zs(sel)), 'b', f(sel)/1e6, imag(Zx(sel)), 'g', 'LineWidth', 1.3); grid on; ylim([-2000 1500]);
xlabel('frequency (MHz)'); ylabel('Im Z (\Omega)');
