%% m16_impedance_kt.m  Plate impedance, fr/fa, and the IEEE Std 176 kt extraction
% Purpose : (1) tabulate fr/fa versus kt from the free-plate resonance condition kt^2 = x cot x,
%           x = (pi/2) fr/fa, and compare the IEEE 176 formula, the small-gap approximation and the
%           lumped keff; (2) "measure" fr and fa on simulated |Z| curves (KLM model) and show how loss
%           and acoustic loading bias the extracted kt.
% Change  : kt list, Qm, tan delta, front/back loads.
% Expected: IEEE 176 returns kt exactly on a lossless free plate (fr/fa = 0.9959, 0.9621, 0.8870,
%           0.7487 for kt = 0.1, 0.3, 0.5, 0.7); the small-gap form is ~6 % high at 0.5; loading biases it.
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
kt176 = @(fr,fa) sqrt(pi/2*fr./fa.*tan(pi/2*(fa-fr)./fa));

%% (1) fr/fa from the resonance condition, and three ways back to k
fprintf('   kt    fr/fa    IEEE176   small-gap   keff\n');
ktl = [0.1 0.3 0.5 0.51 0.6 0.7];
rr  = zeros(size(ktl));
for q = 1:numel(ktl)
    a = 1e-6; b = pi/2 - 1e-12;                           % bisection on x cot x = kt^2
    for it = 1:80, c = (a+b)/2; if c/tan(c) > ktl(q)^2, a = c; else b = c; end, end
    rr(q) = 2*((a+b)/2)/pi;
    fprintf('  %.2f  %.4f   %.4f    %.4f     %.4f\n', ktl(q), rr(q), kt176(rr(q),1), sqrt(pi^2/4*(1-rr(q))), sqrt(1-rr(q)^2));
end

%% (2) simulated impedance curves and extraction
fa = 5e6; base = struct('rho',7500,'v',4600,'kt',0.51,'er',1470,'Qm',Inf,'tand',0,'d',4600/(2*fa),'A',1.5e-6);
f  = (2.5e6:1e3:6.5e6)';
cases = {'lossless, free', Inf, 0, 0, 0; 'Qm 65, tan d 0.02', 65, 0.02, 0, 0; 'Qm 10', 10, 0.02, 0, 0; ...
         'Qm 65, front 1.5 MRayl', 65, 0.02, 0, 1.5; 'Qm 65, back 5 MRayl', 65, 0.02, 5, 0};
fprintf('\n case                       f_m (MHz)  f_n (MHz)  kt(|Z|)   f_s     f_p    kt(G,R)\n');
Zc = zeros(numel(f), size(cases,1));
for q = 1:size(cases,1)
    P = base; P.Qm = cases{q,2}; P.tand = cases{q,3};
    Z = zin(klm(f, P, cases{q,4}*1e6*P.A), cases{q,5}*1e6*P.A); Zc(:,q) = Z;
    [~, im] = min(abs(Z)); [~, iM] = max(abs(Z)); [~, iG] = max(real(1./Z)); [~, iR] = max(real(Z));
    if isinf(P.Qm) && P.tand == 0                         % lossless: Re Z = 0, G and R maxima undefined
        fprintf(' %-26s %8.4f  %8.4f   %6.4f      (no loss: G, R undefined)\n', cases{q,1}, f(im)/1e6, f(iM)/1e6, kt176(f(im),f(iM)));
    else
        fprintf(' %-26s %8.4f  %8.4f   %6.4f  %7.4f %7.4f  %6.4f\n', cases{q,1}, f(im)/1e6, f(iM)/1e6, kt176(f(im),f(iM)), f(iG)/1e6, f(iR)/1e6, kt176(f(iG),f(iR)));
    end
end

%% Figures
figure('Name','impedance curves');
semilogy(f/1e6, abs(Zc(:,1)), 'k', f/1e6, abs(Zc(:,2)), 'b', f/1e6, abs(Zc(:,3)), 'g', f/1e6, abs(Zc(:,4)), 'r', 'LineWidth', 1.2); grid on;
xlabel('frequency (MHz)'); ylabel('|Z| (\Omega)'); title('PZT-5H plate, f_a = 5 MHz, 1.5 mm^2');
legend('lossless free', 'Q_m 65', 'Q_m 10', 'Q_m 65 + 1.5 MRayl front', 'Location', 'southwest');
figure('Name','k from fr/fa');
r = linspace(0.6, 0.9999, 400);
plot(r, kt176(r,1), 'b', r, sqrt(pi^2/4*(1-r)), 'r--', r, sqrt(1-r.^2), 'k:', 'LineWidth', 1.5); grid on;
xlabel('f_r / f_a'); ylabel('coupling'); legend('IEEE 176 k_t', 'small-gap (\pi^2/4)(1 - f_r/f_a)', 'k_{eff} = (1 - (f_r/f_a)^2)^{1/2}');
title('Reading k from the resonance gap');
