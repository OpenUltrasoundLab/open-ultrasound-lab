%% m24_kspace_evanescent.m  --  channel data in (kx, f): the propagating cone, aliasing, evanescent waves
% Purpose : (1) 2D FFT of 0-degree plane-wave channel data. The receive leg is a one-way
%           wave field on the array, so its energy lies inside the cone |kx| <= f/c; the
%           region outside is evanescent and empty. (2) Spatial aliasing: with pitch p the
%           kx axis wraps at 1/(2p); the energy that would fold is measured per target
%           depth for p = 0.30 mm (about lambda at 5 MHz). (3) Why evanescent components are
%           never back-propagated: their gain exp(2 pi z sqrt(kx^2 - (f/c)^2)) explodes.
%           (4) Phase-shift imaging (depth loop, exact plane-wave map, no interpolation)
%           against the Stolt version of the same map (one interpolation): equal images.
% Change  : target depths, pitch, the depth step of the phase-shift loop.
% Expect  : energy outside the cone below -40 dB (window leakage only); the share of
%           energy beyond |kx| = 1/(2 x 0.30 mm) grows from a few % for deep targets to
%           ~30 % for a target at 5 mm; Stolt vs phase-shift residual near -30 dB.
% Model   : as m24_stolt_planewave.m (Born, point elements, no attenuation), pitch 0.15 mm
%           (lambda/2 at 5 MHz, 256 elements) so the reference spectrum is not aliased.
%           Runtime ~10 s in Octave.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);
tStart = tic;

%% Parameters and data (0-degree plane wave, IQ)
c = 1.540; f0 = 5; lam = c/f0; FBW = 0.6; fs = 20;
p = 0.15;  Nel = 256;  xe = ((0:Nel-1) - (Nel-1)/2)*p;        % 38.4 mm aperture
depths = [5 10 20 30];  nt = 1000;
sig_f = FBW*f0/(2*sqrt(2*log(2)));  sig_t = 1/(2*pi*sig_f);
kt = (-ceil(4*sig_t*fs):ceil(4*sig_t*fs)).';  envK = exp(-(kt/fs).^2/(2*sig_t^2));
nfft = 1024;  Df = fs/nfft;  fq = [0:nfft/2-1, -nfft/2:-1]*Df + f0;
bq = find(fq >= 0.8 & fq <= 9.2);  [fB, o] = sort(fq(bq));  bq = bq(o);  fB = fB(:);  nB = numel(fB);
nxp = 512;  kx = [0:nxp/2-1, -nxp/2:-1]/(nxp*p);              % cycles/mm
spec = cell(1, numel(depths));  Sall = 0;
tw = 0.5 - 0.5*cos(pi*min(1, min(0:Nel-1, Nel-1:-1:0)/(0.125*(Nel-1))));   % Tukey(0.25) taper
fprintf('%8s %14s %14s %16s %18s\n', 'depth', 'outside cone', 'outside cone', 'beyond 1/(2x0.3)', 'max receive angle');
fprintf('%8s %14s %14s %16s %18s\n', '(mm)', 'rect (dB)', 'Tukey (dB)', '(% of energy)', 'to array edge');
for id = 1:numel(depths)
    z0 = depths(id);  tau = (z0 + sqrt(xe.^2 + z0^2))/c;      % 0-degree plane wave, two-way
    q = tau*fs;  i0 = floor(q);  fr = q - i0;  w = exp(-1i*2*pi*f0*tau);
    Sd = zeros(nt, Nel);
    Sd(sub2ind([nt Nel], i0 + 1, 1:Nel)) = w.*(1-fr);  Sd(sub2ind([nt Nel], i0 + 2, 1:Nel)) = Sd(sub2ind([nt Nel], i0 + 2, 1:Nel)) + w.*fr;
    S  = conv2(Sd, envK, 'same');  Sall = Sall + S;
    D  = fft(S, nfft, 1);  D = fft(D(bq, :), nxp, 2);  P = abs(D).^2;
    Dw = fft(S, nfft, 1);  Dw = fft(Dw(bq, :).*tw, nxp, 2);  Pw = abs(Dw).^2;   % tapered aperture
    out  = abs(kx) > fB/c;  wrap = abs(kx) > 1/(2*0.30);
    spec{id} = D;
    fprintf('%8.0f %14.1f %14.1f %16.1f %16.1f deg\n', z0, 10*log10(sum(P(out))/sum(P(:))), 10*log10(sum(Pw(out))/sum(Pw(:))), ...
        100*sum(P(wrap & ~out))/sum(P(:)), atan(max(abs(xe))/z0)*180/pi);
end
fprintf('kx Nyquist: %.2f /mm (p = 0.15 mm), %.2f /mm (p = 0.30 mm); cone edge f/c at 5 MHz: %.2f /mm\n', ...
    1/(2*p), 1/(2*0.30), f0/c);
fprintf('the cone reaches the p = 0.30 mm Nyquist at f = c/(2p) = %.2f MHz: above it, steep arrivals alias\n\n', c/(2*0.30));

%% Evanescent gain if one tried to back-propagate
fprintf('back-propagation gain of an evanescent component at 5 MHz, exp(2 pi z sqrt(kx^2-(f/c)^2)):\n');
fprintf('%12s %12s %12s %12s\n', 'kx/(f/c)', 'z = lam/4', 'z = 1 mm', 'z = 10 mm');
for r = [1.01 1.1 1.5]
    a = 2*pi*sqrt(r^2 - 1)*f0/c;                               % decay rate (1/mm)
    fprintf('%12.2f %10.1f dB %10.1f dB %10.0f dB\n', r, 20*log10(exp(a*lam/4)), 20*log10(exp(a*1)), 20*log10(exp(a*10)));
end
fprintf('(the same numbers, negative, are how fast an evanescent field decays away from its source)\n\n');

%% Phase-shift imaging vs Stolt, both with the exact 0-degree map Kz = kz + f/c
D = fft(Sall, nfft, 1);  D = fft(D(bq, :), nxp, 2);            % nB x nxp
dz = 0.077;  zg = (round(3/dz):round(33/dz)).'*dz;  Nz = numel(zg);   % on the Stolt depth grid
KZr = sqrt(max((fB/c).^2 - kx.^2, 0));  prop = (fB/c).^2 >= kx.^2;  % receive kz, propagating mask
t1 = tic;
psi = D .* prop .* exp(1i*2*pi*(KZr + fB/c)*zg(1));            % field at the first depth
step = exp(1i*2*pi*(KZr + fB/c)*dz) .* prop;
imPS = zeros(Nz, nxp);
for k = 1:Nz
    imPS(k, :) = sum(psi, 1);                                  % imaging condition t = 0: sum over f
    psi = psi .* step;                                         % one depth step, transmit + receive
end
imPS = ifft(imPS, [], 2);                                      % kx -> x
tPS = toc(t1);
% Stolt with the same map (interpolation along f instead of the depth loop)
t1 = tic;
ntq = 4096;  Dq = fs/ntq;  fq2 = [0:ntq/2-1, -ntq/2:-1]*Dq + f0;
b2 = find(fq2 >= 0.8 & fq2 <= 9.2);  [fB2, o] = sort(fq2(b2));  b2 = b2(o);  fB2 = fB2(:);  n2 = numel(fB2);
tc = nt/(2*fs);
D2 = fft(Sall, ntq, 1);  D2 = fft(D2(b2, :) .* exp(1i*2*pi*fB2*tc), nxp, 2);
NZ = 2048;  dKz = 1/(NZ*dz);  Kz = (0:NZ-1).'*dKz;  zS = (0:NZ-1).'*dz;
[KX, KZ] = meshgrid(kx, Kz);
F = c*(KX.^2 + KZ.^2)./(2*KZ + eps);  J = abs(c*(KZ.^2 - KX.^2)./(2*KZ.^2 + eps));
ok = KZ > 0 & (KZ - F/c) >= 0;
qq = (F - fB2(1))/Dq + 1;  k0 = floor(qq);  ww = qq - k0;  ok = ok & k0 >= 1 & k0 < n2;
lk = k0(ok) + n2*(ceil(find(ok)/NZ) - 1);
I = zeros(NZ, nxp);  I(ok) = ((1 - ww(ok)).*D2(lk) + ww(ok).*D2(lk+1)) .* J(ok) .* exp(-1i*2*pi*F(ok)*tc);
imS = ifft2(I);  [~, j0] = min(abs(zS - zg(1)));  imS = imS(j0:j0+Nz-1, :);
tS = toc(t1);
% The phase-shift sum over f needs no weight: Stolt's Jacobian df = |df/dKz| dKz is exactly
% what turns the same integral over f into an integral over Kz.
fprintf('phase-shift depth loop: %.2f s (%d depths); Stolt: %.2f s\n', tPS, Nz, tS);
sel = abs(((0:nxp-1)*p + xe(1))) <= 15;  A = imPS(:, sel);  B = imS(:, sel);
al = (A(:)'*B(:))/(A(:)'*A(:));
E1 = abs(al*A);  E2 = abs(B);  m = E2 > 0.01*max(E2(:));
xs = (0:nxp-1)*p + xe(1);  [~, kx0] = min(abs(xs));
fprintf('Stolt (with Jacobian) vs phase-shift (no weight), same pixels:\n');
for z0 = depths
    [~, kz] = min(abs(zg - z0));  a1 = al*imPS(kz, kx0);  b1 = imS(kz, kx0);
    fprintf('  target at %2g mm: amplitude ratio %.3f, phase difference %.3f rad\n', z0, abs(b1)/abs(a1), angle(b1/a1));
end
fprintf('  max |difference| %.1f dB re peak; mean |dB difference| over pixels above -40 dB: %.2f dB\n', ...
    20*log10(max(abs(al*A(:) - B(:)))/max(E2(:))), mean(abs(20*log10(E1(m)./E2(m)))));
fprintf('total run time %.1f s\n', toc(tStart));

%% Figures
figure('Position', [50 50 1100 650]);
for id = [1 3]
    subplot(2,3,1 + (id == 3));
    P = 20*log10(abs(fftshift(spec{id}, 2)) + eps);  P = P - max(P(:));
    imagesc(fftshift(kx), fB, max(P, -50), [-50 0]); axis xy; hold on;
    plot(fB/c, fB, 'r--', -fB/c, fB, 'r--');
    plot([1 1]/(2*0.30), fB([1 end]), 'y:', -[1 1]/(2*0.30), fB([1 end]), 'y:');
    xlabel('k_x (cycles/mm)'); ylabel('f (MHz)'); title(sprintf('|D(k_x,f)| dB, target at %g mm', depths(id)));
end
colormap(gray);
subplot(2,3,3);
kr = linspace(1, 1.5, 200);  plot(kr, -20*log10(exp(2*pi*sqrt(kr.^2 - 1)*f0/c*lam)), 'k'); grid on;
xlabel('k_x / (f/c)'); ylabel('dB per wavelength'); title('evanescent decay over one wavelength, 5 MHz');
xs = (0:nxp-1)*p + xe(1);
subplot(2,3,4); imagesc(xs, zg, max(20*log10(abs(imPS)/max(abs(imPS(:)))), -50), [-50 0]); axis image; axis([-15 15 3 33]);
title('phase-shift (depth loop), dB'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,5); imagesc(xs, zg, max(20*log10(abs(imS)/max(abs(imS(:)))), -50), [-50 0]); axis image; axis([-15 15 3 33]);
title('Stolt (one interpolation), dB'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,6); imagesc(xs, zg, max(20*log10(abs(al*imPS - imS)/max(abs(imS(:)))), -60), [-60 0]); axis image; axis([-15 15 3 33]);
colorbar; title('|difference| (dB re peak)'); xlabel('x (mm)'); ylabel('z (mm)');
