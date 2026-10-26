%% m17_lens_design.m  --  silicone elevation lens: focal length, aberration, apodization, downshift
% Purpose : design a plano-convex RTV-silicone lens for a 40 mm elevation focus,
%           check the paraxial focal length F = R/(c_t/c_l - 1) against an exact
%           Snell ray trace, and quantify what lens attenuation does: a
%           centre-heavy loss (inverse apodization) and a downshift of the
%           pulse-echo spectrum.
% Change  : cl (900-1100 m/s), ct, H (elevation aperture), alpha5, nexp, te.
% Expect  : R = 21.6 mm, sag 0.209 mm, centre 0.509 mm; marginal ray crosses
%           the axis 0.92 mm short of the paraxial focus; two-way loss at 5 MHz
%           3.06 dB (centre) vs 1.80 dB (edge); -6 dB elevation width at focus
%           1.82 mm uniform, 1.80 mm lens-weighted, first sidelobe -26.5 ->
%           -25.7 dB; spectrum centroid 5.00 -> 4.81 MHz; F = 48.0 mm in fat-like
%           1450 m/s tissue.
%           Runtime < 2 s.
clear; close all; rng(1);

%% Parameters
cl = 1000; ct = 1540;          % lens and tissue sound speed (m/s)
Ftarget = 40e-3;               % elevation focus (m)
H = 6e-3;                      % elevation aperture (m)
te = 0.3e-3;                   % lens thickness at the elevation edge (m)
f0 = 5e6;                      % centre frequency (Hz)
alpha5 = 3.0; nexp = 1.5;      % lens attenuation: alpha5 dB/mm at 5 MHz, alpha ~ f^nexp
lam = ct/f0;

%% Geometry
R = Ftarget*(ct/cl - 1);                     % radius of curvature for the target focus
sag = R - sqrt(R^2 - (H/2)^2); Tc = te + sag; % sagitta and centre thickness
T = @(y) Tc - (R - sqrt(R^2 - y.^2));        % lens thickness at elevation y
fprintf('R = %.2f mm, sag = %.3f mm, centre thickness = %.3f mm, F# = %.2f\n', R*1e3, sag*1e3, Tc*1e3, Ftarget/H);

%% Exact ray trace (rays leave the face normally, refract at the spherical surface)
ys = linspace(0.05e-3, H/2, 40); zc = zeros(size(ys));
for k = 1:numel(ys)
  y = ys(k); zs = T(y);                      % surface point (z measured from the face)
  nrm = [y, zs - (Tc - R)]/R;                % outward unit normal (y, z), centre at z = Tc - R
  thi = asin(abs(nrm(1)));                   % incidence angle (ray direction is +z)
  tht = asin(ct/cl*sin(thi));                % Snell: sin(tht)/ct = sin(thi)/cl
  dev = tht - thi;                           % ray turns toward the axis by dev
  zc(k) = zs + y/tan(dev);                   % axis crossing (from the face)
end
Fpar = Tc + R/(ct/cl - 1);                   % paraxial focus measured from the face
fprintf('Paraxial focus %.2f mm from the face (F = %.2f mm from the apex)\n', Fpar*1e3, (Fpar - Tc)*1e3);
fprintf('Marginal ray (y = %.1f mm) crosses at %.2f mm: spherical aberration %.2f mm\n', H/2*1e3, zc(end)*1e3, (zc(end) - Fpar)*1e3);
for ctt = [1450 1540 1580]
  fprintf('  tissue %4d m/s -> F = %.1f mm\n', ctt, R/(ctt/cl - 1)*1e3);
end

%% Attenuation weighting and elevation beam at the focus
alpha = @(fr) alpha5*(fr/5e6).^nexp;                 % dB/mm
y = linspace(-H/2, H/2, 601); dy = y(2) - y(1);
w = 10.^(-alpha(f0)*T(y)*1e3/20);                    % one-way amplitude weight
fprintf('Two-way loss at %.0f MHz: centre %.2f dB, edge %.2f dB\n', f0/1e6, 2*alpha(f0)*Tc*1e3, 2*alpha(f0)*te*1e3);
x = linspace(-4e-3, 4e-3, 801);
Bu = abs(exp(-1i*2*pi*x(:)*y/(lam*Ftarget))*ones(numel(y),1))*dy;
Bw = abs(exp(-1i*2*pi*x(:)*y/(lam*Ftarget))*w(:))*dy;
P2u = (Bu/max(Bu)).^2; P2w = (Bw/max(Bw)).^2;       % two-way (transmit = receive)
wid = @(p) 2*interp1(p(x >= 0), x(x >= 0), 0.5);     % -6 dB full width (two-way amplitude)
fprintf('Two-way -6 dB elevation width: uniform %.2f mm (0.886 lam F/H = %.2f), lens-weighted %.2f mm\n', ...
        wid(P2u)*1e3, 0.886*lam*Ftarget/H*1e3, wid(P2w)*1e3);
pk = @(p) max(p(x > lam*Ftarget/H*0.9 & x < 2.2*lam*Ftarget/H));
fprintf('First sidelobe (two-way): uniform %.1f dB, lens-weighted %.1f dB\n', 20*log10(pk(P2u)), 20*log10(pk(P2w)));

%% Spectral downshift (Gaussian pulse-echo spectrum, 70 % -6 dB bandwidth)
fr = linspace(0.5e6, 10e6, 2000); sig = 0.70*f0/(2*sqrt(2*log(2)));
S0 = exp(-(fr - f0).^2/(2*sig^2));
Tbar = mean(T(y));                                   % aperture-mean thickness
S1 = S0.*10.^(-2*alpha(fr)*Tbar*1e3/20);             % two-way through the mean thickness
cen = @(S) sum(fr.*S)/sum(S);
[~, i1] = max(S1);
fprintf('Mean thickness %.3f mm; centroid %.2f -> %.2f MHz; peak %.2f -> %.2f MHz\n', Tbar*1e3, cen(S0)/1e6, cen(S1)/1e6, f0/1e6, fr(i1)/1e6);
b = 2*alpha5/5*Tbar*1e3/8.686;                        % Np/MHz if alpha were linear (nexp = 1)
fprintf('Linear-law estimate beta*sigma^2 = %.2f MHz (beta = %.3f Np/MHz)\n', b*(sig/1e6)^2, b);

%% Figures
figure('Position', [60 60 1000 620]);
subplot(2,2,1); hold on;
for k = 1:4:numel(ys)
  plot([0 T(ys(k))]*1e3, [ys(k) ys(k)]*1e3, 'b'); plot([T(ys(k)) zc(k)]*1e3, [ys(k) 0]*1e3, 'b');
  plot([0 T(ys(k))]*1e3, -[ys(k) ys(k)]*1e3, 'b'); plot([T(ys(k)) zc(k)]*1e3, -[ys(k) 0]*1e3, 'b');
end
plot(Fpar*1e3, 0, 'ro', 'MarkerFaceColor', 'r'); grid on; xlabel('depth z (mm)'); ylabel('elevation y (mm)');
title(sprintf('Exact rays, R = %.1f mm, c_l = %d m/s', R*1e3, cl)); xlim([0 45]);
subplot(2,2,2); plot(y*1e3, T(y)*1e3, 'LineWidth', 1.5); hold on; plot(y*1e3, w, 'LineWidth', 1.5); grid on;
xlabel('elevation y (mm)'); legend('thickness (mm)', 'one-way weight at f_0', 'Location', 'south'); title('Lens profile and attenuation weight');
subplot(2,2,3); plot(x*1e3, 20*log10(P2u + 1e-6), x*1e3, 20*log10(P2w + 1e-6), 'LineWidth', 1.3); grid on; ylim([-60 0]);
xlabel('elevation x at focus (mm)'); ylabel('two-way (dB)'); legend('uniform', 'lens-weighted'); title('Elevation beam at 40 mm');
subplot(2,2,4); plot(fr/1e6, S0, fr/1e6, S1/max(S1), 'LineWidth', 1.3); grid on; xlabel('frequency (MHz)'); ylabel('normalised');
legend('before lens', 'after lens (two-way, normalised)'); title('Lens loss is a low-pass filter');
