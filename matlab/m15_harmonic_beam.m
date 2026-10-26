%% m15_harmonic_beam.m
% Purpose : why tissue-harmonic beams are narrower and cleaner. Quasi-linear focal-plane model:
%           the second-harmonic source is proportional to the square of the fundamental
%           transmit pattern, and it is received at 2 f0 by the same array.
%             fundamental two-way : [E_f0 AF_f0]_tx * [E_f0 AF_f0]_rx
%             harmonic two-way    : [E_f0 AF_f0]_tx^2 * [E_2f0 AF_2f0]_rx
%             linear at 2 f0      : [E_2f0 AF_2f0]^2   (same probe driven at 2 f0)
% Change  : f0, number of elements, pitch, steering angle, apodization.
% Expect  : continuous apertures (u = x D / (lambda0 F)): -6 dB full widths
%           sinc 1.207, sinc^2 0.886 (x 0.734), sinc^2(u) sinc(2u) 0.511, sinc^2(2u) 0.443;
%           first sidelobes -13.3, -26.5, -27.9, -26.5 dB.
%           Array, 3 MHz, 64 elements, pitch 0.3 mm, steered 30 deg: the 2 f0 receive pattern has a
%           grating lobe near -20.8 deg; in the harmonic two-way pattern it is suppressed by the
%           low transmit level there (see printed levels).
% Runtime : under 1 s.
clear; close all; rng(1);

%% Continuous aperture: widths and sidelobes
u = linspace(0, 6, 600001);
sincf = @(x) (abs(x) < 1e-12) + (abs(x) >= 1e-12).*sin(pi*x)./(pi*x + (abs(x) < 1e-12));
pats = {sincf(u), sincf(u).^2, sincf(u).^2.*sincf(2*u), sincf(2*u).^2};
lab  = {'fundamental one-way sinc', 'fundamental two-way / harmonic tx', 'harmonic two-way', 'linear at 2f0 two-way'};
fprintf('%-36s  -6 dB FW (lambda0 F/D)   peak sidelobe (dB)\n', 'pattern');
for i = 1:4
  a = abs(pats{i}); i6 = find(a < 0.5, 1); d = diff(a); i0 = find(d(1:end-1) < 0 & d(2:end) >= 0, 1) + 1;
  fprintf('%-36s       %6.3f                %7.2f\n', lab{i}, 2*u(i6), 20*log10(max(a(i0:end))));
end

%% Linear array in the focal plane (far-field array factor in sin(theta))
c = 1540; f0 = 3e6; Nel = 64; pitch = 0.3e-3; kerf = 0.03e-3; th0 = 30; apod = 'rect';
w = pitch - kerf; xn = ((0:Nel-1) - (Nel-1)/2)*pitch;
if strcmp(apod, 'hann'), wn = 0.5 - 0.5*cos(2*pi*((0:Nel-1) + 1)/(Nel + 1)); else, wn = ones(1, Nel); end
th = linspace(-90, 90, 3601); s = sind(th); s0 = sind(th0);
AF = @(f) abs(exp(1i*2*pi*f/c*(s(:) - s0)*xn)*wn(:)).'/sum(wn);
EF = @(f) abs(sincf(w*f/c*s));                          % rigid-baffle element factor
P1 = EF(f0).*AF(f0); P2 = EF(2*f0).*AF(2*f0);
fund = P1.^2; harm = P1.^2.*P2; lin2 = P2.^2;
ipk = find(abs(th - th0) < 2);                          % normalize to the steered main lobe
dBn = @(p) 20*log10(p/max(p(ipk)) + 1e-12);
thg1 = asind(max(-1, min(1, s0 - c/f0/pitch))); thg2 = asind(s0 - c/(2*f0)/pitch);
fprintf('\nArray: %d el, pitch %.2f mm, steer %d deg, f0 %.1f MHz\n', Nel, pitch*1e3, th0, f0/1e6);
if abs(s0 - c/f0/pitch) <= 1, g1 = sprintf('%.1f deg', thg1); else, g1 = 'none (|sin| > 1)'; end
fprintf(' predicted grating lobe: tx/rx at f0 %s ; rx at 2f0 %.1f deg\n', g1, thg2);
win = abs(th - thg2) < 3;
names = {'fundamental two-way', 'harmonic two-way', 'linear 2f0 two-way', 'rx only at 2f0'}; P = {fund, harm, lin2, P2};
fprintf(' %-22s  -6 dB width (deg)   level at the 2f0 grating lobe (dB)\n', 'pattern');
for i = 1:4
  d = dBn(P{i}); [~, k0] = min(abs(th - th0)); a = k0; b = k0;   % contiguous -6 dB region around th0
  while a > 1 && d(a-1) >= -6, a = a - 1; end
  while b < numel(th) && d(b+1) >= -6, b = b + 1; end
  fprintf(' %-22s      %6.2f             %7.1f\n', names{i}, th(b) - th(a), max(d(win)));
end

figure('Position', [60 60 900 420]);
plot(th, dBn(fund), 'b', th, dBn(harm), 'r', th, dBn(lin2), 'k--', 'LineWidth', 1.3);
ylim([-80 0]); xlim([-90 90]); grid on; xlabel('angle (deg)'); ylabel('two-way level (dB)');
legend('fundamental', 'harmonic (tx^2 x rx at 2f_0)', 'linear at 2f_0', 'Location', 'southwest');
title(sprintf('Two-way beam, %d elements, pitch %.2f mm, steered %d deg', Nel, pitch*1e3, th0));
