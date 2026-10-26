%% m30_viscoelastic_dispersion.m
% Purpose : Shear-wave dispersion in Kelvin-Voigt and power-law (spring-pot)
%           media: phase velocity, group velocity (numerical d(omega)/dk),
%           attenuation, the dispersion slope (200-600 Hz) and the shape
%           change of a pulse propagated as a 1D plane wave,
%           U(x,w) = U(0,w) exp(-j k(w) x).
% Change  : c0, etaList, n (power-law exponent), fc (pulse centre frequency).
% Expect  : KV, c0 = 1.5 m/s, eta = 1 Pa s: c_ph = 1.54, 1.66, 2.21, 3.16 m/s
%           at 100, 200, 500, 1000 Hz; alpha(200 Hz) = 197 Np/m;
%           c_g(200 Hz) = 2.01 m/s. Low-frequency approximation
%           c0 (1 + 3x^2/8) within 0.1 % at x = 0.3. Power law n = 0.2:
%           c_g/c_ph = 1/(1 - n/2) = 1.111 and alpha = (w/c_ph) tan(n pi/4).
rng(1);
rho = 1000; c0 = 1.5; mu = rho*c0^2;
kKV = @(f, eta) 2*pi*f.*sqrt(rho./(mu + 1j*2*pi*f*eta));
cphKV = @(f, eta) 2*pi*f./real(kKV(f, eta));
h = 0.5; cgKV = @(f, eta) 2*pi*2*h./(real(kKV(f + h, eta)) - real(kKV(f - h, eta)));
fl = [100 200 500 1000];
for eta = [0.4 1 2]
  fprintf('KV eta = %.1f Pa s: c_ph =', eta); fprintf(' %.3f', cphKV(fl, eta));
  fprintf(' m/s; c_g = '); fprintf(' %.3f', cgKV(fl, eta)); fprintf(' m/s; alpha = '); fprintf(' %.0f', -imag(kKV(fl, eta)));
  fprintf(' Np/m; slope(200-600 Hz) = %.2f m/s/kHz\n', (cphKV(600, eta) - cphKV(200, eta))/0.4);
end
x = 0.3; fprintf('Low-frequency check x = 0.3: exact %.4f, 1 + 3x^2/8 = %.4f\n', sqrt(2*(1 + x^2)/(1 + sqrt(1 + x^2))), 1 + 3*x^2/8);
fprintf('High-frequency KV limit sqrt(2 w eta/rho) at 1000 Hz, eta = 1: %.2f m/s (exact %.2f)\n', sqrt(2*2*pi*1000*1/rho), cphKV(1000, 1));
% power law with c_ph(100 Hz) = c0
n = 0.2; w0 = 2*pi*100; cr = c0*cos(n*pi/4);
kPL = @(f) (2*pi*f/cr).*(2*pi*f/w0).^(-n/2).*(cos(n*pi/4) - 1j*sin(n*pi/4));
cphPL = @(f) 2*pi*f./real(kPL(f)); cgPL = @(f) 2*pi*2*h./(real(kPL(f + h)) - real(kPL(f - h)));
fprintf('Power law n = %.1f: c_ph = ', n); fprintf(' %.3f', cphPL(fl)); fprintf(' m/s; c_g/c_ph = '); fprintf(' %.4f', cgPL(fl)./cphPL(fl));
fprintf(' (theory %.4f)\n', 1/(1 - n/2));
fprintf('   alpha(200 Hz) = %.1f Np/m, (w/c_ph) tan(n pi/4) = %.1f Np/m; factor per decade 10^(n/2) = %.3f\n', -imag(kPL(200)), 2*pi*200/cphPL(200)*tan(n*pi/4), 10^(n/2));
% pulse propagation (KV eta = 1, fc = 250 Hz)
fs = 20e3; Nt = 1024; t = (0:Nt-1)'/fs; fc = 250; sf = 0.8*fc/2/sqrt(2*log(2)); st = 1/(2*pi*sf);
u0 = exp(-(t - 0.004).^2/(2*st^2)).*cos(2*pi*fc*(t - 0.004)); U0 = fft(u0);
fq = (0:Nt-1)'*fs/Nt; fq(fq > fs/2) = fq(fq > fs/2) - fs;
xs = [5 10 15]*1e-3; tp = zeros(1,3); U = zeros(Nt,3);
for q = 1:3
  mp = (1:Nt/2-1)'; kp = kKV(mp*fs/Nt, 1); Zp = U0(mp+1).*exp(-1j*kp*xs(q));
  Zf = zeros(Nt,1); Zf(mp+1) = Zp; Zf(Nt-mp+1) = conj(Zp); U(:,q) = real(ifft(Zf));
  Z = fft(U(:,q)); Z(fq < 0) = 0; Z(2:Nt/2) = 2*Z(2:Nt/2); e = abs(ifft(Z));
  [~, im] = max(e); de = e(im-1) - 2*e(im) + e(im+1); tp(q) = (im - 1 + 0.5*(e(im-1) - e(im+1))/de)/fs;
end
fprintf('Pulse fc = %d Hz, KV eta = 1: envelope TTP 5 -> 15 mm gives %.3f m/s (c_ph(fc) = %.3f, c_g(fc) = %.3f)\n', fc, 10e-3/(tp(3) - tp(1)), cphKV(fc, 1), cgKV(fc, 1));
ff = 10:5:1000;
figure(1); subplot(1,2,1); plot(ff, cphKV(ff, 1), 'b', ff, cgKV(ff, 1), 'r--', ff, cphPL(ff), 'k', ff, cgPL(ff), 'k--');
xlabel('f (Hz)'); ylabel('speed (m/s)'); legend('KV c_{ph}', 'KV c_g', 'PL c_{ph}', 'PL c_g', 'location', 'northwest'); title('c_0 = 1.5 m/s, \eta = 1 Pa s, n = 0.2');
subplot(1,2,2); plot(t*1e3, U); xlabel('t (ms)'); ylabel('u'); legend('5 mm', '10 mm', '15 mm'); title('KV pulse, 250 Hz');
