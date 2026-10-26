%% m06_svd_cyst.m  --  Module 06: low-rank approximation of a B-mode cyst phantom
% Purpose : Build a synthetic B-mode image (random scatterers * 2-D RF PSF,
%           anechoic cyst, point targets, layers), detect the envelope with an
%           FFT-based analytic signal, log-compress, then study rank-k
%           approximations: images at ranks 1, 2, 5, 20, 60, the singular-value
%           spectrum, the Eckart-Young identities ||B - B_k||_F^2 = sum_{i>k} sigma_i^2
%           and ||B - B_k||_2 = sigma_{k+1}, a rival rank-k approximation that
%           must lose, and cyst contrast / CNR vs k.
% Change  : f0, PSF widths, cyst radius, dynamic range DR, the list of ranks ks.
% Expect  : both Eckart-Young identities hold (Frobenius to ~1e-14, spectral norm to
%           the printed digits); the rival rank-k approximation loses by 1.5-2.2x;
%           rank 1 already holds 97.3 % of the energy (the -30 dB mean level), so
%           judge the images, not the percentage. Cyst contrast about -21 dB at full
%           rank, capped by the 60 dB display floor (the cyst interior is clipped),
%           and only -6 dB at k = 1; CNR peaks at k = 2 (9.4 vs 3.8 at full rank)
%           because low rank averages speckle, at the price of separable streaks.
% Runtime : ~4 s in Octave.
clear; close all; rng(1);

%% 1. Geometry and RF point-spread function
c0 = 1540; f0 = 5e6; fs = 40e6;              % m/s, Hz, Hz
dz = c0/(2*fs);                              % axial sample spacing (m), 19.25 um
dx = 0.15e-3;                                % lateral line spacing (m)
z = (10e-3:dz:40e-3)'; x = (-19.2e-3:dx:19.2e-3-dx);  % 1559 x 256 grid
Nz = numel(z); Nx = numel(x);
sz = 0.12e-3;                                % axial envelope std (m), FWHM 0.28 mm
sx = 0.30e-3;                                % lateral beam std (m), FWHM 0.71 mm
kz = (-ceil(4*sz/dz):ceil(4*sz/dz))'*dz;  kx = (-ceil(4*sx/dx):ceil(4*sx/dx))*dx;
psf = exp(-kz.^2/(2*sz^2)) .* cos(2*pi*f0*2*kz/c0) * exp(-kx.^2/(2*sx^2));  % RF PSF (two-way phase)

%% 2. Scatterer field: echogenicity map x Gaussian random amplitudes
[X, Z] = meshgrid(x, z);
echo = ones(Nz, Nx);
echo(Z > 14e-3 & Z < 15.5e-3) = 6;           % bright layer (+7.8 dB)
echo(Z > 31e-3 & Z < 33e-3)   = 0.3;         % dark layer (-5.2 dB)
cyst = (X - 2e-3).^2 + (Z - 24e-3).^2 < (4e-3)^2;
echo(cyst) = 0;                              % anechoic cyst, radius 4 mm
incl = ((X + 9e-3)/3e-3).^2 + ((Z - 24e-3)/2.5e-3).^2 < 1;
echo(incl) = 3;                              % hyperechoic inclusion (+4.8 dB)
scat = sqrt(echo) .* randn(Nz, Nx);          % one random scatterer per sample
pts = [-12e-3 20e-3; -4e-3 20e-3; 4e-3 20e-3; 12e-3 20e-3; 10e-3 36e-3];
for p = 1:size(pts,1)
    [~, iz] = min(abs(z - pts(p,2))); [~, ix] = min(abs(x - pts(p,1)));
    scat(iz, ix) = scat(iz, ix) + 400;        % point targets (wires)
end
% 2-D convolution via zero-padded FFTs (base MATLAB conv2 would also work)
Pz = Nz + numel(kz) - 1; Px = Nx + numel(kx) - 1;
rf = real(ifft2(fft2(scat, Pz, Px) .* fft2(psf, Pz, Px)));
oz = (numel(kz)-1)/2; ox = (numel(kx)-1)/2;
rf = rf(oz+(1:Nz), ox+(1:Nx));
rf = rf + 0.002*std(rf(:))*randn(Nz, Nx);    % electronic noise, -54 dB re mean speckle

%% 3. Envelope (analytic signal via FFT, no toolbox hilbert) and log compression
N = Nz; h = zeros(N,1); h(1) = 1;
if mod(N,2) == 0, h(N/2+1) = 1; h(2:N/2) = 2; else, h(2:(N+1)/2) = 2; end
env = abs(ifft(fft(rf) .* repmat(h, 1, Nx)));
env = env(1:4:end, :); zd = z(1:4:end);      % decimate axially to 77 um pixels
cystd = cyst(1:4:end, :);
bgd = ~cystd & abs(Z(1:4:end,:) - 24e-3) < 4e-3 & abs(X(1:4:end,:) - 2e-3) > 6e-3 & abs(X(1:4:end,:) + 9e-3) > 5e-3;
DR = 60;                                     % dynamic range (dB)
B = 20*log10(env/max(env(:)));  B = max(B, -DR);   % log-compressed image (dB), A in the text
[m, n] = size(B);
fprintf('B-mode image %d x %d pixels (%.0f um axial x %.0f um lateral), DR = %d dB\n', m, n, 4*dz*1e6, dx*1e6, DR);

%% 4. SVD, spectrum, Eckart-Young check, cyst metrics vs rank
[U, S, V] = svd(B, 'econ'); s = diag(S);
inC = cystd & ((X(1:4:end,:) - 2e-3).^2 + (Z(1:4:end,:) - 24e-3).^2 < (3e-3)^2);  % cyst core
ks = [1 2 5 10 20 40 80 n];
fprintf('\nEckart-Young: ||B-Bk||_F^2 = sum_{i>k} sigma_i^2 and ||B-Bk||_2 = sigma_{k+1}\n');
fprintf('   k   energy%%   ||B-Bk||_F^2    sum_{i>k} s^2   rel.diff   ||B-Bk||_2  sigma_k+1  contrast(dB)  CNR\n');
con = zeros(size(ks)); cnr = con;
for j = 1:numel(ks)
    k = ks(j);
    Bk = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    e2 = norm(B - Bk, 'fro')^2; tail = sum(s(k+1:end).^2);
    if k < n, e2s = norm(B - Bk, 2); sk1 = s(k+1); else, e2s = 0; sk1 = 0; end
    mc = mean(Bk(inC)); mb = mean(Bk(bgd));
    con(j) = mc - mb;                        % contrast on the dB image
    cnr(j) = abs(mc - mb)/sqrt(var(Bk(inC)) + var(Bk(bgd)));
    if tail > 0, rd = abs(e2 - tail)/tail; else, rd = 0; end
    fprintf('%4d   %7.3f   %13.6e   %13.6e   %8.1e   %9.3f  %9.3f   %7.1f     %5.2f\n', k, 100*sum(s(1:k).^2)/sum(s.^2), e2, tail, rd, e2s, sk1, con(j), cnr(j));
end
fprintf('sigma_1 / sigma_n = %.3g;  ranks for 99%% energy: %d;  for 99.9%%: %d\n', s(1)/s(end), ...
    find(cumsum(s.^2)/sum(s.^2) >= 0.99, 1), find(cumsum(s.^2)/sum(s.^2) >= 0.999, 1));
fprintf('(Energy is dominated by the mean level: judge the residual image, not the percentage.)\n');
% Optimality: any other rank-k matrix does worse. Rival: project B onto the span of
% its k highest-energy columns (column-subset selection, a common "cheap" low rank).
[~, ord] = sort(sum(B.^2, 1), 'descend');
fprintf('\nRival rank-k approximation (span of the k strongest columns) vs truncated SVD, Frobenius error ratio:\n');
for k = [1 5 20 60]
    Q = orth(B(:, ord(1:k))); Bc = Q*(Q'*B);
    fprintf('  k = %2d: ||B-Bc||_F / ||B-Bk||_F = %.3f  (>= 1 by Eckart-Young)\n', k, norm(B - Bc, 'fro')/sqrt(sum(s(k+1:end).^2)));
end

%% 5. Figures
figure('Position', [50 50 1000 640]);
show = [n 1 2 5 20 60]; ttl = {'original (full rank)', 'rank 1', 'rank 2', 'rank 5', 'rank 20', 'rank 60'};
for j = 1:numel(show)
    k = show(j); Bk = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    subplot(2,3,j); imagesc(x*1e3, zd*1e3, Bk, [-DR 0]); colormap(gray); axis image;
    title(ttl{j}); xlabel('lateral (mm)'); if mod(j,3) == 1, ylabel('depth (mm)'); end
end
figure('Position', [50 50 1000 360]);
subplot(1,2,1); semilogy(s/s(1), 'k.-'); grid on; hold on;
semilogy(ks, s(ks)/s(1), 'ro');
xlabel('index i'); ylabel('\sigma_i / \sigma_1'); title('Singular values of the log-compressed image');
subplot(1,2,2); [ax, h1, h2] = plotyy(ks, con, ks, cnr); set(ax, 'XScale', 'log');
set(h1, 'Marker', 'o'); set(h2, 'Marker', 's'); grid on;
xlabel('rank k'); ylabel(ax(1), 'cyst contrast (dB)'); ylabel(ax(2), 'CNR');
title('Cyst metrics of the rank-k image');
