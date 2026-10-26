%% m04_fft_2d_kspace.m  -- Module 04: 2D FFT and the k-space of a pulse-echo image
% Purpose : simulate a small RF image (speckle + point targets) with a broadband,
%           focal-plane (Fraunhofer) pulse-echo PSF, then look at log|FFT2|: the axial band
%           sits at kz = +-2 f0/c, the lateral extent is set by the transmit-receive
%           aperture autoconvolution, |kx| <= f/(c F#) at each frequency (a bow-tie).
%           Also: 2D FFT of a single PSF, of its envelope (B-mode PSF), and of separable windows.
% Change  : f0, FBW, Fnum, apod ('rect' or 'hann'), dx, dz, the scatterer density.
% Expect  : 2f0/c = 6.49 cycles/mm at 5 MHz; two-way lateral edge 1/(lambda F#) = 1.62
%           cycles/mm at F# 2; >99 % of the image k-space energy inside the bow-tie;
%           PSF -6 dB widths close to 0.886 lambda F# (rect) and c tau/2;
%           separable rect window: -13.3 dB sidelobes on the axes, -26.5 dB on the diagonal;
%           filtering: a k-space mask matched to the PSF passband keeps ~99 % of the signal energy
%           and raises the SNR of white noise by ~10 log10(1/mask area); a one-sided kz mask gives
%           the analytic image and its envelope (B-mode).
% Model   : paraxial, focal plane, Tx = Rx aperture, first-order Born, shift-invariant
%           over the small field (circular convolution). Base MATLAB only (R2019b+, Octave 8).
clear; close all; rng(1);

%% Parameters
c    = 1.540;      % speed of sound (mm/us)
f0   = 5;          % center frequency (MHz)
FBW  = 0.60;       % -6 dB fractional bandwidth of the two-way pulse
Fnum = 2;          % F-number z/D (same aperture on transmit and receive)
apod = 'rect';     % 'rect' or 'hann'
Nz = 128; Nx = 128;          % image size (axial x lateral)
dz = 0.02; dx = 0.06;        % pixel size (mm): axial Nyquist 25 c/mm, lateral 8.33 c/mm

B   = FBW*f0; sig_f = B/(2*sqrt(2*log(2)));          % amplitude-spectrum std (MHz)
lam0 = c/f0;
sincf = @(u) (sin(pi*u) + (u == 0)) ./ (pi*u + (u == 0));
if strcmp(apod, 'hann')
    A1 = @(u) sincf(u) + 0.5*sincf(u - 1) + 0.5*sincf(u + 1);   % FT of 0.5+0.5cos(2 pi xi/D), peak 1
else
    A1 = @(u) sincf(u);                                         % FT of a uniform aperture
end

%% Axes (fft order: origin at index 1) and shifted axes for display
kz_f = [0:Nz/2-1, -Nz/2:-1]/(Nz*dz);       % axial spatial frequency (cycles/mm)
x_f  = [0:Nx/2-1, -Nx/2:-1]*dx;            % lateral position (mm), origin first
kz_s = (-Nz/2:Nz/2-1)/(Nz*dz); kx_s = (-Nx/2:Nx/2-1)/(Nx*dx);
z_s  = (-Nz/2:Nz/2-1)*dz;      x_s  = (-Nx/2:Nx/2-1)*dx;

%% PSF: for every temporal frequency f = c kz/2, lateral two-way pattern A1(x/(lambda F#))^2
Hz = zeros(Nz, Nx);                         % axial-FFT x lateral-space representation
for k = 2:Nz/2                              % positive kz (skip DC)
    f  = kz_f(k)*c/2;                       % MHz
    Pf = exp(-(f - f0)^2/(2*sig_f^2));      % two-way pulse amplitude spectrum (on axis)
    u  = x_f/((c/f)*Fnum);                  % normalized lateral coordinate x/(lambda F#)
    Hz(k, :) = Pf * A1(u).^2;               % transmit pattern x receive pattern
end
Hz(Nz:-1:Nz/2+2, :) = Hz(2:Nz/2, :);        % Hermitian twin (real, zero-phase PSF)
h  = real(ifft(Hz, [], 1));                 % RF PSF h(z, x), centered at the origin
ha = ifft([Hz(1, :); 2*Hz(2:Nz/2, :); zeros(Nz/2, Nx)], [], 1);   % analytic PSF (axial)
H  = fft2(h);                               % 2D spectrum of the PSF

%% Scatterers: speckle + three point targets; image = PSF (*) scatterers (circular)
S = 0.5*randn(Nz, Nx) .* (rand(Nz, Nx) < 0.2);  % ~20 scatterers per resolution cell
pts = [32 30; 64 64; 100 100];              % (iz, ix) of the point targets
for ip = 1:size(pts, 1), S(pts(ip, 1), pts(ip, 2)) = 10; end
K   = H .* fft2(S);                         % k-space of the image (2D convolution theorem)
img = real(ifft2(K));
S1 = zeros(Nz, Nx); S1(Nz/2+1, Nx/2+1) = 1; % single point at the center
psf = real(ifft2(H .* fft2(S1)));

%% Measurements
envc = abs(ha); envc = envc/max(envc(:));
ex = fftshift(envc(1, :)); ez = fftshift(envc(:, 1)).';
cross = @(y, xa, lv, i) xa(i-1) + (lv - y(i-1)).*(xa(i) - xa(i-1))./(y(i) - y(i-1));   % first crossing below lv
xr = x_s(Nx/2+1:end); er = ex(Nx/2+1:end); zr = z_s(Nz/2+1:end); fr = ez(Nz/2+1:end);
lat6 = 2*cross(er, xr, 0.5, find(er < 0.5, 1));
ax6  = 2*cross(fr, zr, 0.5, find(fr < 0.5, 1));
tau = 4*log(2)/(pi*B);
kz0 = 2*f0/c; fL = f0 - B/2; fH = f0 + B/2;
kedge = 1/(lam0*Fnum);
[~, ir] = min(abs(kz_f - kz0));                                  % row nearest 2 f0/c
row = abs(H(ir, :)); row = fftshift(row/max(row)); rdB = 20*log10(row + 1e-300);
pos = Nx/2+1:Nx;
k30 = cross(rdB(pos), kx_s(pos), -30, find(rdB(pos) < -30, 1));  % -30 dB half-width
[KX, KZ] = meshgrid(kx_s, kz_s);
inside = abs(KX) <= abs(KZ)/(2*Fnum) + 1/(Nx*dx);                % bow-tie (+1 bin)
E = abs(fftshift(K)).^2;
fprintf('--- k-space of a pulse-echo image (f0 = %.1f MHz, FBW = %.0f %%, F# = %.1f, %s) ---\n', f0, 100*FBW, Fnum, apod);
fprintf('axial band center 2 f0/c = %.3f cycles/mm (RF period lambda/2 = %.4f mm)\n', kz0, lam0/2);
fprintf('axial -6 dB band 2fL/c .. 2fH/c = %.3f .. %.3f cycles/mm\n', 2*fL/c, 2*fH/c);
fprintf('lateral edge at f0: two-way 1/(lambda F#) = %.3f c/mm, one-way %.3f c/mm\n', kedge, kedge/2);
fprintf('measured |H| -30 dB lateral half-width at kz = %.2f c/mm: %.3f c/mm (rect triangle: %.3f)\n', kz_f(ir), k30, kedge*(kz_f(ir)/kz0)*(1 - 10^(-1.5)));
fprintf('fraction of image k-space energy inside the bow-tie |kx| <= |kz|/(2F#): %.5f\n', sum(E(inside))/sum(E(:)));
fprintf('PSF lateral -6 dB width: %.3f mm (narrowband theory %.3f mm = %.3f lambda F#)\n', lat6, ...
        lam0*Fnum*(0.886*strcmp(apod, 'rect') + 1.441*strcmp(apod, 'hann')), 0.886*strcmp(apod, 'rect') + 1.441*strcmp(apod, 'hann'));
fprintf('PSF axial   -6 dB width: %.3f mm (c tau/2 = %.3f mm)\n', ax6, c*tau/2);
fprintf('lateral pixel %.3f mm vs Nyquist for the two-way edge at fH: %.3f mm\n', dx, 1/(2*fH/(c*Fnum)));

%% Figure 1: image, k-space; PSF, k-space
toDB = @(A) max(-60, 20*log10(abs(A)/max(abs(A(:))) + 1e-300));
kzmax = kz_s(end);
figure('Name', 'k-space', 'Position', [60 60 900 720]);
subplot(2, 2, 1); imagesc(x_s, z_s, img); colormap(gray); axis image;
xlabel('lateral x (mm)'); ylabel('axial z (mm)'); title('RF image: speckle + 3 points');
subplot(2, 2, 2); imagesc(kx_s, kz_s, toDB(fftshift(K))); set(gca, 'YDir', 'normal'); hold on;
plot([kx_s(1) kx_s(end)], [kz0 kz0], 'r--', [kx_s(1) kx_s(end)], -[kz0 kz0], 'r--');
plot([0 kzmax/(2*Fnum)], [0 kzmax], 'y--', [0 -kzmax/(2*Fnum)], [0 kzmax], 'y--');
plot([0 kzmax/(2*Fnum)], [0 -kzmax], 'y--', [0 -kzmax/(2*Fnum)], [0 -kzmax], 'y--'); hold off;
xlabel('k_x (cycles/mm)'); ylabel('k_z (cycles/mm)'); title('log|FFT2| (60 dB): red 2f_0/c, yellow bow-tie'); colorbar;
subplot(2, 2, 3); imagesc(x_s, z_s, psf); axis image;
xlabel('lateral x (mm)'); ylabel('axial z (mm)'); title('single RF PSF');
subplot(2, 2, 4); imagesc(kx_s, kz_s, toDB(fftshift(H))); set(gca, 'YDir', 'normal'); hold on;
plot([kx_s(1) kx_s(end)], [kz0 kz0], 'r--', [kx_s(1) kx_s(end)], -[kz0 kz0], 'r--');
plot([-kedge kedge], [kz0 kz0], 'go', 'MarkerFaceColor', 'g'); hold off;
xlabel('k_x (cycles/mm)'); ylabel('k_z (cycles/mm)'); title('PSF k-space; green: \pm1/(\lambda_0 F#)'); colorbar;

%% Separable windows and the B-mode (envelope) PSF
Lw = 16; nw = 0:Lw-1;
wins = {ones(1, Lw), 0.5 - 0.5*cos(2*pi*nw/Lw)};  % rect, periodic Hann
wname = {'rect x rect', 'Hann x Hann'};
fprintf('\n--- separable 2D windows (%d x %d, 8x zero padded) ---\n', Lw, Lw);
figure('Name', 'separable windows and B-mode PSF', 'Position', [60 60 900 720]); colormap(gray);
for iw = 1:2
    w2 = wins{iw}.' * wins{iw};                   % separable: outer product
    W2 = abs(fftshift(fft2(w2, Nz, Nx))); W2 = W2/max(W2(:));
    ax = 20*log10(W2(Nz/2+1, Nx/2+1:end) + 1e-300);              % along kx at kz = 0
    k = 2; while ~(ax(k) < ax(k-1) && ax(k) <= ax(k+1)), k = k + 1; end
    [pslA, iA] = max(ax(k:end)); iA = iA + k - 1;               % first-sidelobe column
    diagv = 20*log10(W2(Nz/2+iA, Nx/2+iA));                      % product sidelobe on the diagonal
    fprintf('%-12s axis sidelobe %.2f dB, diagonal sidelobe %.2f dB (sum of the two 1D levels %.2f dB)\n', ...
            wname{iw}, pslA, diagv, 2*pslA);
    subplot(2, 2, iw); imagesc((-Nx/2:Nx/2-1)/Nx, (-Nz/2:Nz/2-1)/Nz, max(-80, 20*log10(W2 + 1e-300)));
    set(gca, 'YDir', 'normal'); axis image; xlim([-0.25 0.25]); ylim([-0.25 0.25]); colorbar;
    xlabel('k_x (cycles/sample)'); ylabel('k_z (cycles/sample)'); title(['|FFT2| ' wname{iw} ' (dB)']);
end
envp = abs(ifft2(fft2(abs(ha)) .* fft2(S1)));    % B-mode PSF (envelope), centered
subplot(2, 2, 3); imagesc(x_s, z_s, envp); axis image;
xlabel('lateral x (mm)'); ylabel('axial z (mm)'); title('B-mode PSF (envelope)');
subplot(2, 2, 4); imagesc(kx_s, kz_s, toDB(fftshift(fft2(abs(ha))))); set(gca, 'YDir', 'normal'); colorbar;
xlabel('k_x (cycles/mm)'); ylabel('k_z (cycles/mm)'); title('envelope k-space: band moved to baseband');


%% Filtering in k-space: matched passband mask, then a one-sided kz mask for the envelope
Nimg  = randn(Nz, Nx); Nimg = Nimg*norm(img(:))/norm(Nimg(:));   % scaled to exactly 0 dB SNR
Kn    = fft2(Nimg);
Mpass = abs(H) > 0.01*max(abs(H(:)));        % passband: |H| above -40 dB (fft order)
fill  = mean(Mpass(:));                      % fraction of k-space kept
sigF = real(ifft2(K .* Mpass));             % filtered signal part
noiF = real(ifft2(Kn .* Mpass));            % filtered noise part
snr0 = round(100*10*log10(sum(img(:).^2)/sum(Nimg(:).^2)))/100 + 0;
snr1 = 10*log10(sum(sigF(:).^2)/sum(noiF(:).^2));
fprintf('\n--- k-space filtering (white noise at %.1f dB SNR) ---\n', snr0);
fprintf('passband mask |H| > -40 dB keeps %.1f %% of k-space; signal energy retained %.2f %%\n', 100*fill, 100*sum(sigF(:).^2)/sum(img(:).^2));
fprintf('SNR after mask %.1f dB: gain %.1f dB (prediction 10 log10(1/fill) = %.1f dB)\n', snr1, snr1 - snr0, 10*log10(1/fill));
% one-sided axial mask: keep kz > 0 (x2), drop kz < 0 -> analytic image along z, |.| = envelope
Mhalf = [ones(1, Nx); 2*ones(Nz/2-1, Nx); ones(1, Nx); zeros(Nz/2-1, Nx)];
anaI  = ifft2((K + Kn) .* Mpass .* Mhalf);
bmode = 20*log10(abs(anaI)/max(abs(anaI(:))) + 1e-300);
fprintf('one-sided kz mask: max |real(analytic) - filtered RF| = %.1e (the real part is the RF image)\n', ...
        max(max(abs(real(anaI) - (sigF + noiF)))));

figure('Name', 'k-space filtering', 'Position', [60 60 900 720]); colormap(gray);
subplot(2, 2, 1); imagesc(x_s, z_s, img + Nimg); axis image;
xlabel('lateral x (mm)'); ylabel('axial z (mm)'); title(sprintf('RF image + white noise (%.0f dB SNR)', snr0));
subplot(2, 2, 2); imagesc(kx_s, kz_s, toDB(fftshift(K + Kn))); set(gca, 'YDir', 'normal'); hold on;
contour(kx_s, kz_s, double(fftshift(Mpass)), [0.5 0.5], 'g'); hold off; colorbar;
xlabel('k_x (cycles/mm)'); ylabel('k_z (cycles/mm)'); title('noisy k-space; green: passband mask');
subplot(2, 2, 3); imagesc(x_s, z_s, sigF + noiF); axis image;
xlabel('lateral x (mm)'); ylabel('axial z (mm)'); title(sprintf('after the mask: SNR %.1f dB', snr1));
subplot(2, 2, 4); imagesc(x_s, z_s, max(-40, bmode)); axis image; colorbar;
xlabel('lateral x (mm)'); ylabel('axial z (mm)'); title('envelope via one-sided k_z mask (40 dB)');
