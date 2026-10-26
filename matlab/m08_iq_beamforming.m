%% m08_iq_beamforming.m  -  DAS on RF vs DAS on IQ, with and without phase rotation
% Purpose : simulate plane-wave (0 deg) channel data for one point target on
%           a linear array, demodulate each channel to IQ (mix, FIR low-pass,
%           decimate by 4), and beamform three ways:
%             (a) RF DAS (RF upsampled x4, linear interpolation) + analytic
%                 signal along depth  -> reference envelope
%             (b) IQ DAS: interpolate IQ at tau_n, rotate by exp(+j2*pi*fd*tau_n)
%             (c) IQ DAS without the rotation
%           Prints peak loss of (b) and (c) re (a), compares (c) with the
%           narrowband sum |mean(exp(-j2*pi*fd*tau_n))| and the Fresnel
%           integral, and sweeps depth.
% Change  : Nel (32-64), z0, f0, M (IQ decimation), fd (try fd ~= f0).
% Expect  : 32 el, 5 MHz, z = 20 mm: (b) -0.4 dB re (a) (LPF droop plus
%           linear interpolation of IQ at 10 MHz; -0.1 dB with M = 2);
%           (c) -12.4 dB image peak, -13.3 dB on the target (narrowband
%           prediction -13.4 dB, Fresnel -12.5 dB); the loss oscillates with
%           depth and is still -6.5 dB at 60 mm. Runtime ~20 s in Octave.
clear; close all; rng(1);

%% Parameters
c    = 1.54;                 % speed of sound (mm/us)
f0   = 5;  bw = 0.6;         % MHz, -6 dB fractional bandwidth
fd   = f0;                   % demodulation frequency (MHz)
fs   = 40; M = 4; fsIQ = fs/M;   % RF and IQ sampling rates (MHz)
Nel  = 32; pitch = 0.3;      % elements, pitch (mm)
z0   = 20;                   % target depth (mm)
sigf = bw*f0/(2*sqrt(2*log(2)));  sigt = 1/(2*pi*sigf);   % MHz, us
xel  = ((0:Nel-1) - (Nel-1)/2)*pitch;                      % element x (mm), row

% FIR low-pass for the IQ chain (49-tap Hamming windowed sinc, fc = 4 MHz)
L = 49; fcut = 4; k = (-(L-1)/2:(L-1)/2)';
u = 2*fcut/fs*k; sc = (sin(pi*u) + (u==0))./(pi*u + (u==0));
b = (2*fcut/fs)*sc.*(0.54 - 0.46*cos(2*pi*(0:L-1)'/(L-1))); b = b/sum(b);

% Channel data (columns = elements). Plane-wave transmit: t = (z + r)/c
simChannels = @(zt, xt, Nt) exp(-((0:Nt-1)'/fs - (zt + sqrt(zt^2 + (xel - xt).^2))/c).^2/(2*sigt^2)) ...
              .* cos(2*pi*f0*((0:Nt-1)'/fs - (zt + sqrt(zt^2 + (xel - xt).^2))/c));

%% Beamform: point target at (0, z0)
Nt  = ceil((2*z0 + 12)/c*fs) + 64;  Nt = Nt + mod(Nt, 2*M);
t   = (0:Nt-1)'/fs;
rf  = simChannels(z0, 0, Nt);
iqf = 2*conv2(rf.*exp(-1j*2*pi*fd*t), b, 'same');          % mix + LPF (all channels)
iq  = iqf(1:M:end, :);  tIQ = (0:size(iq,1)-1)'/fsIQ;      % decimate
rfU = interpft(rf, 4*Nt); tU = (0:4*Nt-1)'/(4*fs);         % RF upsampled to 160 MHz

xg = linspace(-3, 3, 121);            % lateral grid (mm)
zg = z0 + linspace(-1, 1, 128)';      % depth grid (mm), 15.7 um step
A = zeros(numel(zg), numel(xg)); B = A; C = A;
for ix = 1:numel(xg)
  tau = (zg + sqrt(zg.^2 + (xg(ix) - xel).^2))/c;          % Nz x Nel delays (us)
  sa = zeros(numel(zg),1); sb = sa; scc = sa;
  for n = 1:Nel
    sa  = sa + interp1(tU, rfU(:,n), tau(:,n), 'linear', 0);
    v   = interp1(tIQ, real(iq(:,n)), tau(:,n), 'linear', 0) + ...
          1j*interp1(tIQ, imag(iq(:,n)), tau(:,n), 'linear', 0);
    sb  = sb + v.*exp(1j*2*pi*fd*tau(:,n));
    scc = scc + v;
  end
  A(:,ix) = sa; B(:,ix) = sb; C(:,ix) = scc;
end
% analytic signal of the RF-beamformed column along depth (zero-padded FFT mask)
Np = 2*numel(zg); h = zeros(Np,1); h(1) = 1; h(2:Np/2) = 2; h(Np/2+1) = 1;
Aa = ifft(fft([A; zeros(Np-numel(zg), numel(xg))]).*repmat(h, 1, numel(xg)));
Aenv = abs(Aa(1:numel(zg), :));
pkA = max(Aenv(:)); pkB = max(abs(B(:))); pkC = max(abs(C(:)));

taun = (z0 + sqrt(z0^2 + xel.^2))/c;
pred = 20*log10(abs(mean(exp(-1j*2*pi*fd*taun))));
lam  = c/f0; psi = pi*max(xel)^2/(lam*z0);
s = linspace(0, 1, 4001); fr = 20*log10(abs(trapz(s, exp(-1j*psi*s.^2))));
fprintf('Nel = %d, pitch %.2f mm, f0 = fd = %.1f MHz, z0 = %d mm, IQ rate %.1f MHz\n', Nel, pitch, f0, z0, fsIQ);
fprintf('  (b) IQ + rotation  peak re (a) = %+.2f dB\n', 20*log10(pkB/pkA));
fprintf('  (c) IQ, no rotation peak re (a) = %+.2f dB\n', 20*log10(pkC/pkA));
fprintf('  narrowband prediction |mean(exp(-j2pi fd tau_n))| = %+.2f dB\n', pred);
fprintf('  Fresnel psi = %.2f rad, |int_0^1 exp(-j psi s^2) ds| = %+.2f dB\n', psi, fr);
% lateral -6 dB width of (a) and (b)
latA = max(Aenv)/pkA; latB = max(abs(B))/pkA; latC = max(abs(C))/pkA;
w6 = @(p) (xg(2)-xg(1))*sum(p >= 0.5*max(p));
fprintf('  lateral -6 dB width: (a) %.2f mm, (b) %.2f mm, (c) %.2f mm (lambda F# = %.2f mm)\n', ...
        w6(latA), w6(latB), w6(latC), lam*z0/(Nel*pitch));

%% Depth sweep of the no-rotation loss (peak at the target only)
zs = 10:5:60; lossC = zeros(size(zs)); lossB = lossC; predS = lossC; frS = lossC;
for i = 1:numel(zs)
  zt = zs(i); Nt2 = ceil((2*zt + 12)/c*fs) + 64; Nt2 = Nt2 + mod(Nt2, 2*M); t2 = (0:Nt2-1)'/fs;
  r2 = simChannels(zt, 0, Nt2);
  q2 = 2*conv2(r2.*exp(-1j*2*pi*fd*t2), b, 'same'); q2 = q2(1:M:end,:); tq = (0:size(q2,1)-1)'/fsIQ;
  zz = zt + linspace(-0.3, 0.3, 61)';  tau = (zz + sqrt(zz.^2 + xel.^2))/c;
  sb = zeros(size(zz)); sc2 = sb; sr = sb;
  for n = 1:Nel
    v = interp1(tq, real(q2(:,n)), tau(:,n), 'linear', 0) + 1j*interp1(tq, imag(q2(:,n)), tau(:,n), 'linear', 0);
    sb = sb + v.*exp(1j*2*pi*fd*tau(:,n)); sc2 = sc2 + v;
  end
  tn = (zt + sqrt(zt^2 + xel.^2))/c;
  lossC(i) = 20*log10(max(abs(sc2))/max(abs(sb)));
  predS(i) = 20*log10(abs(mean(exp(-1j*2*pi*fd*tn))));
  ps = pi*max(xel)^2/(lam*zt); frS(i) = 20*log10(abs(trapz(s, exp(-1j*ps*s.^2))));
end
fprintf('\nDepth sweep (on-axis peak of (c) re (b)):\n z (mm)   ');  fprintf('%7d', zs);
fprintf('\n (c)-(b)  '); fprintf('%7.1f', lossC);
fprintf('\n NB sum   '); fprintf('%7.1f', predS);
fprintf('\n Fresnel  '); fprintf('%7.1f', frS); fprintf('\n');

%% Figures
figure('Position', [80 80 950 600]);
dbI = @(I) max(20*log10(abs(I)/pkA + 1e-12), -40);
subplot(2,3,1); imagesc(xg, zg, dbI(Aenv)); caxis([-40 0]); colormap(gray);
xlabel('x (mm)'); ylabel('z (mm)'); title('(a) RF DAS');
subplot(2,3,2); imagesc(xg, zg, dbI(B)); caxis([-40 0]); xlabel('x (mm)'); title('(b) IQ + rotation');
subplot(2,3,3); imagesc(xg, zg, dbI(C)); caxis([-40 0]); xlabel('x (mm)'); title('(c) IQ, no rotation'); colorbar;
subplot(2,3,4:5);
plot(xg, 20*log10(latA), 'k', xg, 20*log10(latB), 'b--', xg, 20*log10(latC), 'r'); ylim([-50 3]);
xlabel('lateral x (mm)'); ylabel('dB re (a) peak'); legend('(a) RF', '(b) IQ + rot', '(c) IQ no rot');
title('lateral profile (max over depth)'); grid on;
subplot(2,3,6); plot(zs, lossC, 'ro-', zs, predS, 'k.--', zs, frS, 'b:'); grid on;
xlabel('depth (mm)'); ylabel('loss (dB)'); legend('sim', 'NB sum', 'Fresnel', 'Location', 'southeast');
title('no-rotation peak loss vs depth');
