%% m15_attenuation_dispersion.m
% Purpose : power-law attenuation alpha(f) = alpha0 f^y, the dispersion that causality
%           (Kramers-Kronig, nearly local form) attaches to it, the downshift of a
%           Gaussian pulse, and the pulse after propagation (FFT model).
% Change  : alpha0 (dB/cm/MHz), y, depth z, f0, FBW.
% Expect  : y = 1, alpha0 = 0.5: c(10 MHz) - c(1 MHz) = 3.19 m/s (1.38 m/s per neper of ln f)
%           y = 1.1: 3.56 m/s ; y = 1.5: 4.71 m/s
%           5 MHz, 60 % pulse, 0.5 dB/cm/MHz: centroid shift 0.467 MHz per 5 cm one way
%           (sigma_f^2 * alpha0' * z), 0.93 MHz for a 5 cm round trip
% Runtime : under 1 s.
clear; close all; rng(1);

%% Dispersion from Kramers-Kronig (Szabo 1994 / Waters et al. 2000 nearly local form)
c0 = 1540;                 % phase speed at the reference frequency f_ref (m/s)
fref = 1e6;                % reference frequency (Hz)
f = linspace(0.5e6, 15e6, 600);
dBperNp = 20*log10(exp(1));          % 8.686 dB per neper
cKK = @(a0dB, y, fq) 1./( 1/c0 + ...
     (abs(y-1) < 1e-6) .* ( -2/pi*(a0dB/dBperNp*100/(2*pi*1e6)^y) .* log(fq/fref) ) + ...
     (abs(y-1) >= 1e-6) .* ( (a0dB/dBperNp*100/(2*pi*1e6)^y) .* tan(pi*y/2) .* ((2*pi*fq).^(y-1) - (2*pi*fref).^(y-1)) ) );
cases = [0.5 1.0; 0.5 1.1; 0.5 1.5; 1.0 1.0; 0.5 1.0001];
fprintf('alpha0 (dB/cm/MHz)   y     c(5)-c(1)   c(10)-c(1)  (m/s)\n');
for i = 1:size(cases,1)
  a0 = cases(i,1); y = cases(i,2);
  fprintf('   %5.2f            %6.4f   %7.2f     %7.2f\n', a0, y, cKK(a0,y,5e6) - c0, cKK(a0,y,10e6) - c0);
end
fprintf('y -> 1 limit check: the y = 1.0001 row must equal the y = 1 row.\n');

%% Gaussian pulse downshift: amplitude spectrum exp(-(f-f0)^2/(2 sf^2)) * exp(-alpha0' f z)
f0 = 5e6; FBW = 0.6; sf = FBW*f0/(2*sqrt(2*log(2)));       % sf = 1.274 MHz
a0 = 0.5; a0p = a0/dBperNp*100/1e6;                          % Np/(m*Hz)
fg = linspace(0, 15e6, 3001);
fprintf('\n5 MHz, 60 %% pulse, alpha0 = 0.5 dB/cm/MHz, y = 1 (one-way):\n');
fprintf(' z (cm)  centroid (MHz)  theory f0 - sf^2 a0 z (MHz)   loss at f0 (dB)\n');
for zc = [0 2.5 5 10]
  z = zc/100; S = exp(-(fg - f0).^2/(2*sf^2)).*exp(-a0p*fg*z);
  fcen = sum(fg.*S)/sum(S);
  fprintf(' %5.1f      %6.3f          %6.3f                        %5.1f\n', zc, fcen/1e6, (f0 - sf^2*a0p*z)/1e6, a0*f0/1e6*zc);
end
fprintf('Round trip to 5 cm (10 cm path): shift %.2f MHz.\n', sf^2*a0p*0.10/1e6);

%% Pulse after propagation (time domain via FFT), y = 1.1 for a visible dispersion
y = 1.1; a0 = 0.5; zc = 8; z = zc/100;
fs = 100e6; N = 4096; t = ((0:N-1) - N/2)/fs;
st = 1/(2*pi*sf); s = exp(-t.^2/(2*st^2)).*cos(2*pi*f0*t);
F = fft(ifftshift(s)); fk = (0:N-1)*fs/N; fk(fk > fs/2) = fk(fk > fs/2) - fs;
af = abs(fk)/1e6; alpha = a0/dBperNp*100*af.^y;              % Np/m
cf = cKK(a0, y, max(abs(fk), 1e3)); cref = cKK(a0, y, f0);
ph = -2*pi*fk*z.*(1./cf - 1/cref);                             % phase relative to a frame moving at c(f0)
G = F.*exp(-alpha*z).*exp(1i*ph);
r = fftshift(real(ifft(G)));
% centroid and envelope peak
R = abs(fft(ifftshift(r))); half = fk > 0;
fprintf('\nAfter %d cm with y = %.1f: spectral centroid %.2f MHz, peak amplitude %.3f of input\n', zc, y, ...
  sum(fk(half).*R(half).^1)/sum(R(half))/1e6, max(abs(r)));
figure('Position', [60 60 900 600]);
subplot(2,2,1); plot(f/1e6, cKK(0.5,1,f) - c0, 'b', f/1e6, cKK(0.5,1.1,f) - c0, 'r', f/1e6, cKK(1,1,f) - c0, 'k', 'LineWidth', 1.3);
xlabel('f (MHz)'); ylabel('c(f) - c(1 MHz) (m/s)'); grid on; legend('0.5, y=1', '0.5, y=1.1', '1.0, y=1', 'Location', 'southeast');
title('Kramers-Kronig dispersion');
subplot(2,2,2); hold on;
for zc2 = [0 5 10]
  S = exp(-(fg - f0).^2/(2*sf^2)).*exp(-a0p*fg*zc2/100); plot(fg/1e6, S, 'LineWidth', 1.3);
end
xlabel('f (MHz)'); ylabel('amplitude'); legend('0 cm', '5 cm', '10 cm'); grid on; title('Downshift, 0.5 dB/cm/MHz');
subplot(2,1,2); plot(t*1e6, s, 'Color', [0.6 0.6 0.6]); hold on; plot(t*1e6, r/max(abs(r)), 'b', 'LineWidth', 1.3);
xlim([-1.2 1.2]); xlabel('t - z/c(f_0) (\mus)'); ylabel('normalized p');
legend('input', sprintf('after %d cm (normalized)', zc)); grid on; title('Attenuated, dispersed pulse');
