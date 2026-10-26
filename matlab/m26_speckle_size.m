%% m26_speckle_size.m  --  speckle size from the PSF: autocovariance widths vs f0, bandwidth, F#
% Purpose : build fully developed speckle IQ images in k-space (white circular Gaussian
%           scatterer field times the pulse-echo transfer function H), then measure the
%           FWHM of the intensity autocovariance axially and laterally and compare with
%             axial   (Gaussian pulse, B = -6 dB two-way bandwidth): 0.8825/B in time,
%                     i.e. 0.441 c/B in depth, equal to the pulse-envelope FWHM
%             lateral (focal plane, CW, rect tx and rx apertures, two-way PSF sinc^2):
%                     0.835 lambda F#   (two-way -6 dB beam width is 0.886 lambda F#)
%           Theory = FWHM of |mu|^2, where mu = inverse FT of |H|^2 (Siegert relation).
%           Also prints the envelope and log-image ACF widths, which are narrower.
% Change  : cases (f0, fractional bandwidth, F#, apodization), grid sizes.
% Expect  : 5 MHz, 60 %, F#2, rect: axial 0.227 mm, lateral 0.514 mm (theory); measured
%           within a few percent; log-image axial width about 0.77 x the intensity width.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);
c  = 1540;                                     % sound speed (m/s)
Nz = 256;  Nx = 256;  dz = 20e-6;  dx = 40e-6; % grid: 5.1 mm x 10.2 mm
uz = [0:Nz/2-1, -Nz/2:-1]' / (Nz*dz);          % axial spatial frequency (cycles/m), FFT order
ux = [0:Nx/2-1, -Nx/2:-1]  / (Nx*dx);          % lateral spatial frequency (cycles/m)
nreal = 4;                                     % independent images averaged per case
cases = [5 0.6 2 0;  3 0.6 2 0;  7.5 0.6 2 0;  5 0.3 2 0;  5 0.9 2 0;  5 0.6 1 0;  5 0.6 3 0;  5 0.6 2 1];
% columns: f0 (MHz), fractional -6 dB bandwidth, F#, apodization (0 rect, 1 Hann)

% FWHM of a sampled, peak-normalized symmetric profile p(1) = 1 (index 1 = zero lag)
fwhm = @(p) 2*((find(p < 0.5, 1) - 2) + (p(find(p < 0.5, 1) - 1) - 0.5) / (p(find(p < 0.5, 1) - 1) - p(find(p < 0.5, 1))));
acf2 = @(I) real(ifft2(abs(fft2(I - mean(I(:)))).^2));     % circular autocovariance

% Lateral transfer function: autocorrelation of the apodization (coarray), sampled on ux.
% The aperture coordinate is scaled so that the coarray half-width is 1/(lambda F#).
coarray = @(u, U, apod) interp1(linspace(-1, 1, 2001), ...
    conv(apod(linspace(-0.5, 0.5, 1001)), apod(linspace(-0.5, 0.5, 1001))) / ...
    sum(apod(linspace(-0.5, 0.5, 1001)).^2), u / U, 'linear', 0);
rectw = @(t) ones(size(t));
hannw = @(t) 0.5 + 0.5*cos(2*pi*t);

fprintf(' f0   BW   F#  apod |  axial FWHM (mm) th / meas   | lateral FWHM (mm) th / meas  (th/lamF#) | env/I  log/I (axial) | speckles/mm^2\n');
for k = 1:size(cases, 1)
    f0 = cases(k,1)*1e6;  B = cases(k,2)*f0;  Fn = cases(k,3);  lam = c/f0;
    if cases(k,4) == 0, apod = rectw; an = 'rect'; else, apod = hannw; an = 'hann'; end
    sf = B / (2*sqrt(2*log(2)));               % |H(f)| = exp(-f^2/(2 sf^2)), -6 dB full width = B
    Hz = exp(-(c*uz/2).^2 / (2*sf^2));         % axial: temporal f = c*uz/2 (round trip)
    Hx = coarray(ux, 1/(lam*Fn), apod);        % lateral: coarray, half-width 1/(lambda F#)
    H  = Hz * Hx;                              % separable pulse-echo transfer function
    % theory: |mu|^2 with mu = IFFT(|H|^2), on a 16x finer 1D grid
    Lf = 4096;  uzf = [0:Lf/2-1, -Lf/2:-1]'/(Lf*dz/4);  muz = abs(ifft(exp(-(c*uzf/2).^2/sf^2)));
    muz = (muz/muz(1)).^2;  thz = fwhm(muz) * dz/4;
    uxf = [0:Lf/2-1, -Lf/2:-1]/(Lf*dx/4);  mux = abs(ifft(coarray(uxf, 1/(lam*Fn), apod).^2));
    mux = (mux/mux(1)).^2;  thx = fwhm(mux) * dx/4;
    Vc = 1 / (sum(muz)*dz/4 * sum(mux)*dx/4);  % 1 / speckle correlation area (int int |mu|^2), per m^2
    % simulation
    AI = 0;  AA = 0;  AL = 0;
    for r = 1:nreal
        W = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
        z = ifft2(fft2(W) .* H);  I = abs(z).^2;
        AI = AI + acf2(I);  AA = AA + acf2(sqrt(I));  AL = AL + acf2(log(I));
    end
    pz = AI(:,1)/AI(1,1);  px = AI(1,:)'/AI(1,1);
    mz = fwhm(pz)*dz;  mx = fwhm(px)*dx;
    ea = fwhm(AA(:,1)/AA(1,1))*dz / mz;  la = fwhm(AL(:,1)/AL(1,1))*dz / mz;
    fprintf('%4.1f %4.1f %3.0f  %s |      %6.3f / %6.3f         |      %6.3f / %6.3f   (%5.3f)    |  %4.2f   %4.2f        | %6.2f\n', ...
        f0/1e6, B/f0, Fn, an, thz*1e3, mz*1e3, thx*1e3, mx*1e3, thx/(lam*Fn), ea, la, Vc*1e-6);
    if k == 1, Ishow = I; pz1 = pz; px1 = px; mzT = muz; mxT = mux; end
end
fprintf('axial theory 0.8825/B * c/2 = 0.441 c/B; lateral rect theory 0.835 lambda F#\n');

%% Figure: one speckle image and its autocovariance profiles vs theory
figure('Position', [80 80 1100 380]);
subplot(1,3,1);
imagesc((0:Nx-1)*dx*1e3, (0:Nz-1)*dz*1e3, 10*log10(Ishow/mean(Ishow(:))), [-30 10]);
colormap(gray); axis image; colorbar; xlabel('lateral (mm)'); ylabel('depth (mm)');
title('speckle, 5 MHz, 60 %, F#2 (dB)');
subplot(1,3,2);
lagz = (0:Nz/2)*dz*1e3;  plot(lagz, pz1(1:Nz/2+1), 'bo'); hold on;
lf = (0:Lf/2)*dz/4*1e3;  plot(lf, mzT(1:Lf/2+1), 'b-', 'LineWidth', 1.5);
plot([0 1], [0.5 0.5], 'k--'); xlim([0 0.6]); grid on;
xlabel('axial lag (mm)'); ylabel('intensity autocovariance'); title('axial: sim (o) vs |mu|^2');
subplot(1,3,3);
lagx = (0:Nx/2)*dx*1e3;  plot(lagx, px1(1:Nx/2+1), 'ro'); hold on;
lf = (0:Lf/2)*dx/4*1e3;  plot(lf, mxT(1:Lf/2+1), 'r-', 'LineWidth', 1.5);
plot([0 2], [0.5 0.5], 'k--'); xlim([0 1.5]); grid on;
xlabel('lateral lag (mm)'); ylabel('intensity autocovariance'); title('lateral: sim (o) vs |mu|^2');
