%% m24_stolt_planewave.m  --  Stolt f-k migration for plane waves versus DAS
% Purpose : simulate complex baseband (IQ) channel data of a 128-element linear array for
%           steered plane waves (point targets, and an anechoic cyst in speckle), then form
%           the image three ways on the same pixel grid:
%             DAS   : delay and sum, dynamic receive aperture F# 1.5 (rectangular)
%             ERM   : Stolt f-k migration with the exploding-reflector fit of Garcia et al.
%                     (2013): v = c/sqrt(1+cos(th)+sin(th)^2), depth factor
%                     beta = (1+cos th)^1.5/(1+cos th+sin^2 th), lateral shear sin th/(2-cos th)
%             exact : the same pipeline with the exact transmit + receive wavevector map
%                     f = c(Kx^2+Kz^2) / (2(Kx sin th + Kz cos th))
%           Pipeline (both maps): FFT in t, FFT in x, receive fan filter equivalent to F# 1.5,
%           2x lateral upsampling, steering compensation exp(j2 pi f x sin(th)/c), ONE linear
%           interpolation along f per Kx column, Jacobian |df/dKz|, IFFT2.
%           Printed: time per angle, -6 dB widths, peak sidelobe, peak position error at
%           0 and 9 deg, cyst contrast, clutter beside the speckle block, residual against
%           DAS, and how interpolation, zero-padding and time-centering change the result.
% Change  : thetaDeg (compounding set), Fnum, ntFFT, tCenter, useFan, pitch.
% Expect  : widths equal to DAS within ~2 % (lateral 0.55 mm single angle, 0.50 mm for
%           7 angles; axial 0.24 mm), PSL within ~1.5 dB, cyst contrast within ~1 dB.
%           ERM at 9 deg: peak displaced laterally by tens of um; exact map: < 5 um.
%           Linear interpolation on a 2048-point FFT (no centering) smears every strong
%           echo laterally (clutter beside the block -11 dB instead of -19 dB); 4096 points
%           with time-centering removes it. f-k is a few times faster than DAS in Octave.
% Model   : linear acoustics, first-order Born, infinite plane-wave transmit, point-like
%           elements (no directivity), no attenuation or spreading, Gaussian two-way pulse
%           (60 % band). Runtime ~20 s in Octave.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);
tStart = tic;

%% Parameters
c     = 1.540;                 % speed of sound (mm/us)
f0    = 5;                     % center frequency (MHz)
lam   = c/f0;                  % wavelength (mm) = 0.308
FBW   = 0.60;                  % -6 dB fractional bandwidth (amplitude spectrum)
fs    = 20;                    % IQ sampling rate (MHz)
Nel   = 128; pitch = 0.30;     % elements, pitch (mm)
xe    = ((0:Nel-1) - (Nel-1)/2) * pitch;   % element positions (mm)
Fnum  = 1.5;                   % receive F-number (DAS aperture z/Fnum, f-k fan filter)
alphaMax = atan(1/(2*Fnum));   % receive half-angle allowed by the F-number (rad)
useFan = true;                 % receive fan filter in f-k (the k-space F-number)
thetaDeg = -9:3:9;             % plane-wave angles (deg); must contain 0 and 9
theta = thetaDeg*pi/180;  Nang = numel(theta);  iMid = find(thetaDeg == 0);
density = 40;                  % speckle scatterers per mm^2
nt    = 1250;                  % samples per channel (62.5 us), t = 0 at transmit
ntFFT = 4096;                  % time FFT length for f-k (3.3x zero-padding)
tCenter = nt/(2*fs);           % time-centering reference (us) before f interpolation

sig_f = FBW*f0/(2*sqrt(2*log(2)));          % Gaussian spectrum std (MHz) = 1.27
sig_t = 1/(2*pi*sig_f);                     % envelope std (us) = 0.125
kt    = (-ceil(4*sig_t*fs):ceil(4*sig_t*fs)).';
envK  = exp(-(kt/fs).^2/(2*sig_t^2));       % sampled envelope

%% Phantoms: (A) point targets, (B) speckle with an anechoic cyst
ptsX = [-10 0 10 0].';  ptsZ = [15 15 15 22].';  ptsA = ones(4,1);
R = 3; cx = 0; cz = 33;  xr = [-8 8]; zr = [27 39];
Ns  = round(density*diff(xr)*diff(zr));
spX = xr(1) + diff(xr)*rand(Ns,1);  spZ = zr(1) + diff(zr)*rand(Ns,1);
keep = (spX-cx).^2 + (spZ-cz).^2 > R^2;  spX = spX(keep);  spZ = spZ(keep);
spA = (randn(numel(spX),1) + 1i*randn(numel(spX),1))/sqrt(2);
fprintf('speckle scatterers %d, cyst radius %.0f mm, %d angles (%s deg)\n', numel(spX), R, Nang, mat2str(thetaDeg));

%% Image grids (fixed): lateral 2x upsampled, depth from a fixed Kz grid
nxp  = 256;  U = 2;  nxf = U*nxp;  dxf = pitch/U;
xf   = xe(1) + (0:nxf-1)*dxf;               % image x grid (mm), 0.15 mm
kxd  = [0:nxp/2-1, -nxp/2:-1]/(nxp*pitch);  % data wavenumbers (cycles/mm)
Kx   = [0:nxf/2-1, -nxf/2:-1]/(nxf*dxf);    % image wavenumbers (cycles/mm)
NZ   = 1024;  dKz = 2*(fs/2048)/c;          % Kz step (cycles/mm): depth extent 78.8 mm
Kz   = (0:NZ-1).'*dKz;  zf = (0:NZ-1).'/(NZ*dKz);   % image z grid (mm), 0.077 mm
[KX, KZ] = meshgrid(Kx, Kz);
ix = find(abs(xf) <= 12);  iz = find(zf >= 10 & zf <= 40);
xg = xf(ix);  zg = zf(iz);  [XG, ZG] = meshgrid(xg, zg);
fprintf('pixel grid %d x %d (dx %.3f mm, dz %.4f mm)\n', numel(zg), numel(xg), dxf, zg(2)-zg(1));

%% DAS geometry (angle independent)
Px = XG(:);  Pz = ZG(:);  nPix = numel(Px);
Rrx = sqrt((Px - xe).^2 + Pz.^2)/c;          % receive time (us), nPix x Nel
act = abs(xe - Px) <= Pz/(2*Fnum);           % dynamic aperture
[pp, ee] = find(act);  rxT = Rrx(act);  clear Rrx
Wn  = 1 ./ accumarray(pp, 1, [nPix 1]);      % unit-sum rectangular weights
Wc  = Wn(pp) .* exp(1i*2*pi*f0*rxT);         % weight x receive part of the remodulation
txTime = @(x, z, th) (x*sin(th) + z*cos(th))/c;

% f-k jobs: [map (1 ERM, 2 exact), FFT length, nearest (1) or linear (0), time-centering]
jobMain = [1 ntFFT 0 1; 2 ntFFT 0 1];
jobTest = [2 2048 1 0; 2 2048 0 0; 2 4096 0 0; 2 4096 0 1; 2 8192 0 1];

%% Simulate, then beamform every angle three ways
img  = zeros(numel(zg), numel(xg), Nang, 3, 2);     % (z, x, angle, method, phantom)
imgT = zeros(numel(zg), numel(xg), size(jobTest,1), 2);
tDAS = 0;  tFK = zeros(1,2);
for ia = 1:Nang
    th = theta(ia);
    S = cell(1,2);
    for ph = 1:2                                 % channel data (IQ, nt x Nel)
        if ph == 1, X0 = ptsX; Z0 = ptsZ; A0 = ptsA; else, X0 = spX; Z0 = spZ; A0 = spA; end
        tau = txTime(X0, Z0, th) + sqrt((X0 - xe).^2 + Z0.^2)/c;   % two-way delay (us)
        w   = A0 .* exp(-1i*2*pi*f0*tau);        % exact carrier phase
        q   = tau*fs;  i0 = floor(q);  fr = q - i0;
        li  = i0 + 1 + nt*(0:Nel-1);
        Sd  = accumarray(li(:), w(:).*(1-fr(:)), [nt*Nel 1]) + accumarray(li(:)+1, w(:).*fr(:), [nt*Nel 1]);
        S{ph} = conv2(reshape(Sd, nt, Nel), envK, 'same');
    end
    % ---- DAS
    t1 = tic;
    T  = txTime(Px(pp), Pz(pp), th) + rxT;  q = T*fs;  i0 = floor(q);  fr = q - i0;
    li = i0 + 1 + nt*(ee - 1);  txRot = exp(1i*2*pi*f0*txTime(Px, Pz, th));
    for ph = 1:2
        v = ((1-fr).*S{ph}(li) + fr.*S{ph}(li+1)) .* Wc;
        img(:,:,ia,1,ph) = reshape(accumarray(pp, v, [nPix 1]) .* txRot, numel(zg), numel(xg));
    end
    tDAS = tDAS + toc(t1);
    % ---- f-k
    jobs = jobMain;  if ia == iMid, jobs = [jobMain; jobTest]; end
    for jb = 1:size(jobs, 1)
        mode = jobs(jb,1);  ntq = jobs(jb,2);  nearest = jobs(jb,3);  tc = jobs(jb,4)*tCenter;
        t1 = tic;
        Dfq  = fs/ntq;  fq = [0:ntq/2-1, -ntq/2:-1]*Dfq + f0;     % RF frequency of each bin
        bq   = find(fq >= 0.8 & fq <= 9.2);  [fB, o] = sort(fq(bq));  bq = bq(o);
        fB   = fB(:);  nB = numel(fB);
        if mode == 1        % Garcia ERM fit
            v  = c/sqrt(1 + cos(th) + sin(th)^2);  be = (1 + cos(th))^1.5/(1 + cos(th) + sin(th)^2);
            ga = sin(th)/(2 - cos(th));  kzh = (KZ - ga*KX)/be;
            F  = v*sqrt(KX.^2 + kzh.^2);  J = v*kzh./(be*sqrt(KX.^2 + kzh.^2) + eps);  ok = kzh > 0;
        else                % exact transmit + receive wavevector map
            M  = KX*sin(th) + KZ*cos(th);  N2 = KX.^2 + KZ.^2;
            F  = c*N2./(2*M + eps);  J = abs(c*(2*KZ.*M - N2*cos(th))./(2*M.^2 + eps));
            ok = M > 0 & (KZ - F*cos(th)/c) >= 0;
        end
        q  = (F - fB(1))/Dfq + 1;
        if nearest, i0 = round(q); w = zeros(size(q)); else, i0 = floor(q); w = q - i0; end
        ok = ok & i0 >= 1 & i0 < nB;
        li = i0(ok) + nB*(ceil(find(ok)/NZ) - 1);  w = w(ok);
        Wt = J(ok) .* exp(-1i*2*pi*F(ok)*tc);    % Jacobian, undo the time-centering
        for ph = 1:2
            D  = fft(S{ph}, ntq, 1);
            D  = D(bq, :) .* exp(1i*2*pi*fB*tc);         % RF frequency rows, time-centered
            D  = fft(D, nxp, 2);                         % data wavenumber kx_d
            if useFan, D = D .* (abs(kxd) <= fB*sin(alphaMax)/c); end
            Du = zeros(nB, nxf);  Du(:, 1:nxp/2) = D(:, 1:nxp/2);  Du(:, end-nxp/2+1:end) = D(:, nxp/2+1:end);
            d  = ifft(Du, [], 2) .* exp(1i*2*pi*(fB*sin(th)/c)*xf);   % steering compensation
            D2 = fft(d, [], 2);                          % nB x nxf, image Kx
            I  = zeros(NZ, nxf);
            I(ok) = ((1 - w).*D2(li) + w.*D2(li+1)) .* Wt;
            im = ifft2(I);
            if jb <= 2, img(:,:,ia,mode+1,ph) = im(iz, ix); else, imgT(:,:,jb-2,ph) = im(iz, ix); end
        end
        if jb <= 2, tFK(mode) = tFK(mode) + toc(t1); end
    end
end
fprintf('time per angle (2 phantoms): DAS %.2f s, f-k ERM %.2f s, f-k exact %.2f s -> f-k %.1fx faster\n', ...
    tDAS/Nang, tFK(1)/Nang, tFK(2)/Nang, tDAS/mean(tFK));
fprintf('work per angle and phantom: DAS %d (pixel, element) pairs; f-k FFTs %dx%d + 2 x %dx%d, %d interpolations\n\n', ...
    numel(pp), ntFFT, Nel, nB, nxf, nnz(ok));

%% Normalize each method by its peak on the central point target (0 deg)
names = {'DAS', 'f-k ERM', 'f-k exact'};
nrm = @(I3) max(max(abs(I3(abs(zg-15) < 1, abs(xg) < 1))));
for m = 1:3, sc = nrm(img(:,:,iMid,m,1));  img(:,:,:,m,:) = img(:,:,:,m,:)/sc; end
for k = 1:size(jobTest,1), sc = nrm(imgT(:,:,k,1));  imgT(:,:,k,:) = imgT(:,:,k,:)/sc; end

%% Metrics on the central target (0, 15 mm)
xl = -1.5:0.005:1.5;  zl = -1.2:0.002:1.2;
box = abs(XG) <= 4 & abs(ZG - 15) <= 4;  mlobe = abs(XG) <= 1.0 & abs(ZG - 15) <= 0.45;
sets = {iMid, 1:Nang, find(thetaDeg == 9)};  setName = {'0 deg', sprintf('%d angles', Nang), '9 deg'};
fprintf('%-10s %-9s %8s %8s %7s %9s %9s %8s\n', 'method', 'set', 'lat-6dB', 'ax-6dB', 'PSL', 'dx peak', 'dz peak', 'A(22)/');
fprintf('%-10s %-9s %8s %8s %7s %9s %9s %8s\n', '', '', '(mm)', '(mm)', '(dB)', '(um)', '(um)', 'A(15)');
for s = 1:3
    for m = 1:3
        E  = abs(sum(img(:,:,sets{s},m,1), 3));
        sub = E .* box;  [pk, k] = max(sub(:));  [kz, kx] = ind2sub(size(E), k);
        rz = max(1,kz-12):min(numel(zg),kz+12);  rx = max(1,kx-8):min(numel(xg),kx+8);
        Lc = interp1(zg(rz), E(rz, kx), zg(kz) + zl, 'spline');  [~, b] = max(Lc);  zpk = zg(kz) + zl(b);
        rowv = interp1(zg(rz), E(rz, rx), zpk, 'spline');
        Lr = interp1(xg(rx), rowv, xg(kx) + xl, 'spline');  [pr, a] = max(Lr);  xpk = xg(kx) + xl(a);
        colv = interp1(xg(rx), E(rz, rx).', xpk, 'spline');
        Lc = interp1(zg(rz), colv(:), zg(kz) + zl, 'spline');  [pc, b] = max(Lc);  zpk = zg(kz) + zl(b);
        wl = sum(Lr >= pr/2)*0.005;  wa = sum(Lc >= pc/2)*0.002;
        PSL = 20*log10(max(E(box & ~mlobe))/pk);
        A22 = max(max(E(abs(ZG-22) < 1 & abs(XG) < 1)))/pk;
        fprintf('%-10s %-9s %8.3f %8.3f %7.1f %9.1f %9.1f %8.2f\n', names{m}, setName{s}, wl, wa, PSL, xpk*1e3, (zpk-15)*1e3, A22);
    end
end
fprintf('theory: receive-only 1.207 lam F# = %.3f mm; axial -6 dB of the envelope 2.355 sig_t c/2 = %.3f mm\n', ...
    1.207*lam*Fnum, 2*sqrt(2*log(2))*sig_t*c/2);
fprintf('pixel grid: dx %.0f um, dz %.0f um (peak positions are spline-refined)\n\n', dxf*1e3, (zg(2)-zg(1))*1e3);

%% Speckle: cyst contrast, clutter beside the block, residual against DAS
rr   = sqrt((XG-cx).^2 + (ZG-cz).^2);
cin  = rr <= R - 0.5;  cout = rr >= R + 1 & rr <= R + 3.5 & ZG <= 38.5;
blk  = abs(XG) < 7 & ZG > 28 & ZG < 38 & rr > R + 1;   side = XG > 9.5 & ZG > 29 & ZG < 38;
gres = @(I, ref) 20*log10(norm((sum(conj(I).*ref, 2)./sum(abs(I).^2, 2)).*I - ref, 'fro')/norm(ref, 'fro'));
spk  = @(I) [10*log10(mean(abs(I(cin)).^2)/mean(abs(I(cout)).^2)), 10*log10(mean(abs(I(side)).^2)/mean(abs(I(blk)).^2))];
fprintf('%-10s %-9s %10s %10s %12s\n', 'method', 'set', 'contrast', 'beside', 'residual vs');
fprintf('%-10s %-9s %10s %10s %12s\n', '', '', '(dB)', 'block (dB)', 'DAS (dB)');
for s = 1:2
    ref = sum(img(:,:,sets{s},1,1) + img(:,:,sets{s},1,2), 3);
    for m = 1:3
        Ic = sum(img(:,:,sets{s},m,1) + img(:,:,sets{s},m,2), 3);  r = spk(sum(img(:,:,sets{s},m,2), 3));
        if m == 1, rs = '(reference)'; else, rs = sprintf('%.1f', gres(Ic, ref)); end
        fprintf('%-10s %-9s %10.1f %10.1f %12s\n', names{m}, setName{s}, r(1), r(2), rs);
    end
end
fprintf('(residual after one complex gain per depth row: removes the depth-gain difference\n');
fprintf(' between unit-sum DAS and Stolt; what remains is the k-space weighting difference)\n\n');

%% Interpolation along f, zero-padding and time-centering (0 deg, exact map)
fprintf('%-8s %6s %8s %9s %10s %10s %12s\n', 'interp', 'ntFFT', 'df(kHz)', 'centered', 'contrast', 'beside', 'residual vs');
fprintf('%-8s %6s %8s %9s %10s %10s %12s\n', '', '', '', '', '(dB)', 'block (dB)', 'DAS (dB)');
ref = img(:,:,iMid,1,1) + img(:,:,iMid,1,2);
for k = 1:size(jobTest,1)
    Ic = imgT(:,:,k,1) + imgT(:,:,k,2);  r = spk(imgT(:,:,k,2));
    nm = 'linear'; if jobTest(k,3), nm = 'nearest'; end
    yn = 'no'; if jobTest(k,4), yn = 'yes'; end
    fprintf('%-8s %6d %8.2f %9s %10.1f %10.1f %12.1f\n', nm, jobTest(k,2), fs/jobTest(k,2)*1e3, yn, r(1), r(2), gres(Ic, ref));
end
fprintf('\ntotal run time %.1f s\n', toc(tStart));

%% Figures
dr = 50;  dB = @(I) max(20*log10(abs(I)/max(abs(I(:))) + eps), -dr);
figure('Position', [40 40 1150 720]);
for m = 1:3
    Ic = sum(img(:,:,:,m,1) + img(:,:,:,m,2), 3);
    subplot(2,3,m); imagesc(xg, zg, dB(Ic), [-dr 0]); axis image; colormap(gray);
    title(sprintf('%s, %d angles (dB)', names{m}, Nang)); xlabel('x (mm)'); ylabel('z (mm)');
end
subplot(2,3,4);
E1 = abs(img(:,:,iMid,1,1));  E2 = abs(img(:,:,iMid,3,1));  E3 = abs(sum(img(:,:,:,3,1),3));
[~, jz] = min(abs(zg - 15));
plot(xg, 20*log10(E1(jz,:)/max(E1(jz,:))), 'k', xg, 20*log10(E2(jz,:)/max(E2(jz,:))), 'r--', ...
     xg, 20*log10(E3(jz,:)/max(E3(jz,:))), 'b');
axis([-12 12 -60 0]); grid on; xlabel('x (mm)'); ylabel('dB');
legend('DAS, 0 deg', 'f-k exact, 0 deg', 'f-k exact, 7 angles'); title('lateral profile, z = 15 mm');
subplot(2,3,5);
Ic = imgT(:,:,2,1) + imgT(:,:,2,2);
imagesc(xg, zg, dB(Ic), [-dr 0]); axis image; title('0 deg, linear, 2048 pts, no centering');
xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,6);
Ic = imgT(:,:,4,1) + imgT(:,:,4,2);
imagesc(xg, zg, dB(Ic), [-dr 0]); axis image; title('0 deg, linear, 4096 pts, centered');
xlabel('x (mm)'); ylabel('z (mm)');
