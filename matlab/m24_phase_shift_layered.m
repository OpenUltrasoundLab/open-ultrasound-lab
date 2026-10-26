%% m24_phase_shift_layered.m  --  phase-shift migration through a layer vs Stolt with one speed
% Purpose : a 0-degree plane wave images point targets below an 8 mm fat-like layer
%           (c1 = 1450 m/s) on tissue (c2 = 1540 m/s). Channel data use ray-traced (Fermat)
%           receive delays and the exact vertical transmit delay. Three reconstructions:
%             phase-shift : depth loop with the local speed, step exp(j2 pi (kz(z) + f/c(z)) dz),
%                           kz(z) = sqrt((f/c(z))^2 - kx^2); kx is conserved across a flat
%                           interface (Snell), so this is exact for horizontal layers
%             Stolt 1540  : exact 0-degree map with c = 1540 m/s everywhere
%             Stolt cbar  : the same with the time-average speed down to 25 mm
%           Printed: peak depth error, -6 dB lateral width and peak amplitude per target.
% Change  : c1, layer thickness d1, target depths, Fnum.
% Expect  : phase-shift peaks within a few um of the truth, lateral width 0.555 mm at
%           every depth. Stolt at 1540 m/s: every target 0.50 mm too deep (theory
%           (c2/c1 - 1) d1 = 0.497 mm) and ~15 % wider. Stolt at cbar = 1510 m/s: depth exact
%           at 25 mm only (+0.20 mm at 15 mm, -0.20 mm at 35 mm).
% Model   : Born, point elements, no attenuation, no transmission losses or multiples at
%           the interface, IQ at 20 MHz, pitch 0.30 mm, fan filter F# 1.5.
%           Runtime ~3 s in Octave.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);
tStart = tic;

%% Parameters
c2 = 1.540; c1 = 1.450; d1 = 8;               % mm/us, mm
f0 = 5; FBW = 0.6; fs = 20; lam = c2/f0;
Nel = 128; pitch = 0.30; xe = ((0:Nel-1) - (Nel-1)/2)*pitch;
Fnum = 1.5; alphaMax = atan(1/(2*Fnum));
ptsX = [0 0 0 6].';  ptsZ = [15 25 35 25].';  nt = 1200;
cz = @(z) c1*(z < d1) + c2*(z >= d1);          % speed profile
sig_f = FBW*f0/(2*sqrt(2*log(2)));  sig_t = 1/(2*pi*sig_f);
kt = (-ceil(4*sig_t*fs):ceil(4*sig_t*fs)).';  envK = exp(-(kt/fs).^2/(2*sig_t^2));

%% Channel data: vertical transmit, Fermat receive path through the interface
Ttx = d1/c1 + (ptsZ - d1)/c2;                  % us, exact for a normal plane wave
xa = repmat(xe, numel(ptsX), 1);  xb = repmat(ptsX, 1, Nel);   % bracket for the crossing point
lo = min(xa, xb);  hi = max(xa, xb);  Zt = repmat(ptsZ, 1, Nel);  Xt = xb;  Xe = xa;
for it = 1:60                                  % bisection on dT/dxr (monotonic in xr)
    xr = (lo + hi)/2;
    g  = (xr - Xe)./(c1*sqrt((xr - Xe).^2 + d1^2)) - (Xt - xr)./(c2*sqrt((Xt - xr).^2 + (Zt - d1).^2));
    hi(g > 0) = xr(g > 0);  lo(g <= 0) = xr(g <= 0);
end
Trx = sqrt((xr - Xe).^2 + d1^2)/c1 + sqrt((Xt - xr).^2 + (Zt - d1).^2)/c2;
tau = Ttx + Trx;                               % two-way delay (us), targets x elements
q = tau*fs;  i0 = floor(q);  fr = q - i0;  w = exp(-1i*2*pi*f0*tau);
li = i0 + 1 + nt*(0:Nel-1);
Sd = accumarray(li(:), w(:).*(1-fr(:)), [nt*Nel 1]) + accumarray(li(:)+1, w(:).*fr(:), [nt*Nel 1]);
S  = conv2(reshape(Sd, nt, Nel), envK, 'same');
% Snell check: the ray to the edge element of the 15 mm target
[~, k] = max(abs(xe));
s1 = (xr(1,k) - xe(k))/sqrt((xr(1,k) - xe(k))^2 + d1^2);  s2 = (ptsX(1) - xr(1,k))/sqrt((ptsX(1) - xr(1,k))^2 + (ptsZ(1) - d1)^2);
fprintf('Snell check (15 mm target, edge element): sin(a1)/c1 = %.5f, sin(a2)/c2 = %.5f (1/(mm/us))\n', s1/c1, s2/c2);

%% Common spectra and grids
nfft = 2048;  Df = fs/nfft;  fq = [0:nfft/2-1, -nfft/2:-1]*Df + f0;
bq = find(fq >= 0.8 & fq <= 9.2);  [fB, o] = sort(fq(bq));  bq = bq(o);  fB = fB(:);  nB = numel(fB);
nxp = 256;  kx = [0:nxp/2-1, -nxp/2:-1]/(nxp*pitch);  xs = xe(1) + (0:nxp-1)*pitch;
D = fft(S, nfft, 1);  D = fft(D(bq, :), nxp, 2);
D = D .* (abs(kx) <= fB*sin(alphaMax)/c2);     % receive fan, F# 1.5 in tissue
dz = 0.0385;  zg = (0:round(40/dz)).'*dz;  Nz = numel(zg);

%% Phase-shift migration through the layered profile
t1 = tic;  psi = D;  imPS = zeros(Nz, nxp);
stp = cell(1,2);                               % one depth step in each medium
for m = 1:2
    cl = c1*(m == 1) + c2*(m == 2);
    kzr = sqrt(max((fB/cl).^2 - kx.^2, 0));  prop = (fB/cl).^2 >= kx.^2;
    stp{m} = exp(1i*2*pi*(kzr + fB/cl)*dz) .* prop;   % receive kz + transmit f/c
end
for k = 1:Nz
    imPS(k, :) = sum(psi, 1);                  % imaging condition: sum over f
    psi = psi .* stp{1 + (cz(zg(k) + dz/2) == c2)};
end
imPS = ifft(imPS, [], 2);  tPS = toc(t1);

%% Stolt (exact 0-degree map) with a single speed
cbar = 25/(d1/c1 + (25 - d1)/c2);              % time-average speed to 25 mm
ntq = 4096;  Dq = fs/ntq;  fq2 = [0:ntq/2-1, -ntq/2:-1]*Dq + f0;
b2 = find(fq2 >= 0.8 & fq2 <= 9.2);  [fB2, o] = sort(fq2(b2));  b2 = b2(o);  fB2 = fB2(:);  n2 = numel(fB2);
tc = nt/(2*fs);
D2 = fft(S, ntq, 1);  D2 = fft(D2(b2, :) .* exp(1i*2*pi*fB2*tc), nxp, 2);
NZ = 2048;  Kz = (0:NZ-1).'/(NZ*dz);  [KX, KZ] = meshgrid(kx, Kz);
imS = cell(1,2);  cs = [c2 cbar];
for m = 1:2
    c = cs(m);  Dm = D2 .* (abs(kx) <= fB2*sin(alphaMax)/c);
    F = c*(KX.^2 + KZ.^2)./(2*KZ + eps);  J = abs(c*(KZ.^2 - KX.^2)./(2*KZ.^2 + eps));
    ok = KZ > 0 & (KZ - F/c) >= 0;  qq = (F - fB2(1))/Dq + 1;  k0 = floor(qq);  ww = qq - k0;
    ok = ok & k0 >= 1 & k0 < n2;  lk = k0(ok) + n2*(ceil(find(ok)/NZ) - 1);
    I = zeros(NZ, nxp);  I(ok) = ((1 - ww(ok)).*Dm(lk) + ww(ok).*Dm(lk+1)) .* J(ok) .* exp(-1i*2*pi*F(ok)*tc);
    im = ifft2(I)/(dz*Df);  imS{m} = im(1:Nz, :);   % same scale as the phase-shift sum over f bins
end

%% Metrics
names = {'phase-shift', 'Stolt 1540', sprintf('Stolt %.0f', cbar*1e3)};
ims = {imPS, imS{1}, imS{2}};
fprintf('time-average speed to 25 mm: %.1f m/s; expected depth error at 1540 m/s for a target at z:\n', cbar*1e3);
fprintf('   dz = (c2/c1 - 1) d1 = %.3f mm (one-way time preserved), same for every target below the layer\n\n', (c2/c1 - 1)*d1);
fprintf('%-12s %8s %10s %9s %9s\n', 'method', 'target', 'dz peak', 'lat-6dB', 'peak');
fprintf('%-12s %8s %10s %9s %9s\n', '', '(x,z)', '(um)', '(mm)', '(dB)');
xq = -1.5:0.005:1.5;  refpk = zeros(1, numel(ptsX));
for m = 1:3
    E = abs(ims{m});  if m == 1, nrm = max(E(:)); end;  E = E/nrm;
    for it = 1:numel(ptsX)
        sx = abs(xs - ptsX(it)) <= 3;  sz = abs(zg - ptsZ(it)) <= 2;
        sub = E(sz, sx);  [pk, k] = max(sub(:));  [kz, kx2] = ind2sub(size(sub), k);
        zz = zg(sz);  xx = xs(sx);
        zfine = zz(kz) + (-0.2:0.001:0.2);  Lc = interp1(zz, sub(:, kx2), zfine, 'spline');  [pk, b] = max(Lc);
        rowv = interp1(zz, sub, zfine(b), 'spline');  Lr = interp1(xx, rowv, xx(kx2) + xq, 'spline');
        wl = sum(Lr >= max(Lr)/2)*0.005;
        fprintf('%-12s %8s %10.0f %9.3f %9.1f\n', names{m}, sprintf('%g,%g', ptsX(it), ptsZ(it)), (zfine(b) - ptsZ(it))*1e3, wl, 20*log10(pk));
    end
end
fprintf('\nphase-shift loop: %.2f s for %d depth steps (cost ~ Nz Nf Nkx); a Stolt pass: one interpolation\n', tPS, Nz);
fprintf('total run time %.1f s\n', toc(tStart));

%% Figures
figure('Position', [50 50 1100 420]);
for m = 1:3
    subplot(1,3,m); E = abs(ims{m});
    imagesc(xs, zg, max(20*log10(E/max(E(:))), -50), [-50 0]); axis image; axis([-12 12 5 40]); colormap(gray); hold on;
    plot([-12 12], [d1 d1], 'y--'); plot(ptsX, ptsZ, 'r+');
    title(sprintf('%s (dB); red + = truth', names{m})); xlabel('x (mm)'); ylabel('z (mm)');
end
