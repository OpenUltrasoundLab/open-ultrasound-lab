%% m23_planewave_compounding.m  --  coherent plane-wave compounding, simulated end to end
% Purpose : a self-contained plane-wave simulator. It synthesizes complex baseband (IQ)
%           channel data of a 128-element linear array for (A) five point targets,
%           (B) a speckle phantom with an anechoic cyst and (C) electronic noise, for 31
%           steered plane waves over +-15 deg. Each angle is beamformed with DAS
%           (dynamic receive aperture, F# 1.5); the complex images are then summed over
%           nested subsets of 1, 3, 11 and 31 angles (same +-15 deg span, steps 30, 15, 3
%           and 1 deg). Printed per subset: angle step, transmit-replica distance
%           lambda/d(sin), equivalent transmit F#, -6 dB lateral width, peak sidelobe and
%           integrated clutter around the central target, cyst contrast / CNR / gCNR,
%           speckle SNR0, and the SNR gain against electronic noise. DAS weights are
%           normalized to unit sum per pixel, so the gain does not change with depth.
% Change  : thetaMax_deg and the subsets (try 11 angles over +-5 deg), Fnum, rxApod
%           ('rect' or 'hann'), FBW, density, the noise is unit-variance per channel sample.
% Expect  : SNR gain within 0.1 dB of 10 log10(N) = 0, 4.8, 10.4, 14.9 dB.
%           1 angle: receive-only focusing; -6 dB width 0.57 mm, close to the one-way
%           1.207 lambda F# = 0.56 mm; peak sidelobe -17 dB (the X-shaped edge arms).
%           31 angles: width 0.46 mm (narrowband theory for a +-15 deg span: 0.45 mm),
%           peak sidelobe -30 dB, cyst contrast -16 -> -38 dB, gCNR 0.87 -> 0.995.
%           Transmit replicas sit at lambda/d(sin theta): 1.19 mm for the 15 deg step
%           (a broad -25 dB shoulder), 5.95 mm for the 3 deg step (about -55 dB, because
%           the receive pattern is already far down there) and 17.9 mm for 1 deg.
%           Speckle SNR0 stays ~1.8-1.9 for every N: coherent compounding sharpens the
%           PSF, it does not average speckle.
% Model   : linear acoustics, first-order Born, ideal plane-wave transmit limited to the
%           array's projected footprint (no edge waves), no attenuation or spreading,
%           point-like elements (no directivity), Gaussian two-way pulse. IQ channel data
%           at fs = 20 MHz: each echo gets its exact carrier phase and a linearly
%           interpolated position, then is convolved with the sampled envelope (envelope
%           std broadened by ~1 %). Runtime ~15-25 s in Octave.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);
tic;

%% Parameters
c     = 1.540;                 % speed of sound (mm/us)
f0    = 5;                     % center frequency (MHz)
lam   = c/f0;                  % wavelength (mm) = 0.308
FBW   = 0.60;                  % -6 dB fractional bandwidth of the two-way pulse spectrum
fs    = 20;                    % IQ sampling rate (MHz)
Nel   = 128; pitch = 0.30;     % elements, pitch (mm)
xe    = ((0:Nel-1) - (Nel-1)/2) * pitch;   % element x positions (mm), row vector
L     = Nel*pitch;             % aperture length (mm) = 38.4
Fnum  = 1.5;                   % receive F-number (dynamic aperture z/Fnum)
rxApod = 'rect';               % 'rect' or 'hann' receive apodization
thetaMax_deg = 15;             % plane-wave span is +-thetaMax
Nang  = 31;                    % angles in the full set (1 deg step)
theta = linspace(-thetaMax_deg, thetaMax_deg, Nang) * pi/180;
iMid  = (Nang+1)/2;
subsets = {iMid, [1 iMid Nang], 1:3:Nang, 1:Nang};     % nested: 1, 3, 11, 31 angles
density = 70;                  % speckle scatterers per mm^2 (about 7 per resolution cell)

sig_f = FBW*f0/(2*sqrt(2*log(2)));        % std of the Gaussian amplitude spectrum (MHz)
sig_t = 1/(2*pi*sig_f);                   % std of the Gaussian envelope (us) = 0.125
kt    = (-ceil(4*sig_t*fs):ceil(4*sig_t*fs)).';
envK  = exp(-(kt/fs).^2/(2*sig_t^2));     % sampled envelope (column), peak 1

%% Phantoms
% (A) point targets: a central target for the metrics, four others for context
ptsX = [0 -5 5 -5 5].';  ptsZ = [30 25 25 35 35].';  ptsA = ones(5,1);
% (B) speckle with an anechoic cyst of radius R at (cx, cz)
R = 3; cx = 0; cz = 30;
xr = [-8.6 8.6]; zr = [21.4 38.6];
Ns = round(density*diff(xr)*diff(zr));
spX = xr(1) + diff(xr)*rand(Ns,1);  spZ = zr(1) + diff(zr)*rand(Ns,1);
keep = (spX-cx).^2 + (spZ-cz).^2 > R^2;
spX = spX(keep); spZ = spZ(keep);
spA = (randn(numel(spX),1) + 1i*randn(numel(spX),1))/sqrt(2);
fprintf('speckle scatterers: %d (%.0f per mm^2), cyst radius %.1f mm\n', numel(spX), density, R);

%% Image grids and angle-independent receive geometry
dx = 0.15; dz = 0.10;                    % B-mode grid (mm), above the IQ Nyquist limits
xg = -8:dx:8;  zg = 22:dz:38;
[XG, ZG] = meshgrid(xg, zg);
xl = -8:0.02:8;                          % fine lateral line through the central target
Px = [XG(:); xl(:)];  Pz = [ZG(:); 30*ones(numel(xl),1)];
nImg = numel(XG);  nPix = numel(Px);

Rrx  = sqrt((Px - xe).^2 + Pz.^2) / c;   % nPix x Nel, receive time (us)
uap  = (xe - Px) ./ (Pz/(2*Fnum));       % position in the receive aperture, |u| <= 1
if strcmp(rxApod, 'hann'), W = 0.5 + 0.5*cos(pi*uap); else, W = ones(size(uap)); end
act  = abs(uap) <= 1;                    % active (pixel, element) pairs only
W    = W .* act;  W = W ./ sum(W, 2);    % normalize: same gain at every depth
[pp, ee] = find(act);
rxT  = Rrx(act);                          % receive time of each active pair
Wc   = W(act) .* exp(1i*2*pi*f0*rxT);     % apodization x receive part of the re-modulation
clear Rrx uap W act
fprintf('active (pixel, element) pairs: %d (%.0f per pixel)\n', numel(pp), numel(pp)/nPix);

txTime = @(x, z, th) (x*sin(th) + z*cos(th))/c;          % plane-wave arrival time (us)
inBeam = @(x, z, th) abs(x - z*tan(th)) <= L/2;          % projected footprint of the array
RsA = sqrt((ptsX - xe).^2 + ptsZ.^2)/c;  EA = exp(-1i*2*pi*f0*RsA);   % receive legs, targets
RsB = sqrt((spX  - xe).^2 + spZ.^2)/c;   EB = exp(-1i*2*pi*f0*RsB);   % receive legs, speckle

%% Simulate and beamform every angle
imgA = zeros(nPix, Nang); imgB = zeros(nPix, Nang); imgC = zeros(nPix, Nang);
for ia = 1:Nang
    th  = theta(ia);
    Ttx = txTime(Px, Pz, th);                          % transmit time to each pixel (us)
    T   = Ttx(pp) + rxT;                               % two-way DAS delay of each pair
    t0  = min([min(T), min(RsB(:)) + min(txTime(spX, spZ, th))]) - 1;   % time origin (us)
    nt  = ceil((max([max(T), max(RsB(:)) + max(txTime(spX, spZ, th))]) + 1 - t0)*fs) + 2;
    S = cell(1, 3);
    for ph = 1:2                                       % channel data, IQ, nt x Nel
        if ph == 1, ttx = txTime(ptsX, ptsZ, th); Rs = RsA; E = EA; a = ptsA .* inBeam(ptsX, ptsZ, th);
        else,       ttx = txTime(spX,  spZ,  th); Rs = RsB; E = EB; a = spA  .* inBeam(spX,  spZ,  th); end
        w  = (a .* exp(-1i*2*pi*f0*ttx)) .* E;         % exact baseband phase of every echo
        q  = (ttx + Rs - t0)*fs;  i0 = floor(q);  fr = q - i0;
        li = i0 + (1 + nt*(0:Nel-1));                  % linear index of the earlier sample
        w1 = w.*(1-fr);  w2 = w.*fr;                     % linear split between two samples
        Sd = accumarray(li(:), w1(:), [nt*Nel 1]) + accumarray(li(:)+1, w2(:), [nt*Nel 1]);
        S{ph} = conv2(reshape(Sd, nt, Nel), envK, 'same');   % pulse envelope
    end
    S{3} = conv2((randn(nt, Nel) + 1i*randn(nt, Nel))/sqrt(2), envK, 'same');  % noise, same band
    q  = (T - t0)*fs;  i0 = floor(q);  fr = q - i0;
    li = i0 + 1 + nt*(ee - 1);
    txRot = exp(1i*2*pi*f0*Ttx);                       % transmit part of the re-modulation
    for ph = 1:3                                       % DAS: interpolate, weight, sum
        v = ((1-fr).*S{ph}(li) + fr.*S{ph}(li+1)) .* Wc;
        img = accumarray(pp, v, [nPix 1]) .* txRot;
        if ph == 1, imgA(:,ia) = img; elseif ph == 2, imgB(:,ia) = img; else, imgC(:,ia) = img; end
    end
end
fprintf('simulation + DAS of %d angles: %.1f s\n\n', Nang, toc);

%% Compound and measure
inB   = 1:nImg;  inL = nImg+1:nPix;
dXc   = XG(:);   dZc = ZG(:) - 30;                         % offsets from the central target
box   = abs(dXc) <= 4 & abs(dZc) <= 4;
mlobe = abs(dXc) <= 1.2*lam*Fnum & abs(dZc) <= 0.4;        % mainlobe region (excluded)
rr    = sqrt((XG(:)-cx).^2 + (ZG(:)-cz).^2);
cin   = rr <= R - 0.5;  cout = rr >= R + 1 & rr <= R + 3.5;
nS    = numel(subsets);
res   = zeros(nS, 14);  prof = zeros(numel(xl), nS);  bA = cell(1,nS);  bB = cell(1,nS);
fprintf('%4s %6s %8s %6s %7s %7s %7s %7s %8s %9s %6s %6s %5s %8s %7s\n', 'N', 'step', 'replica', 'ghost', 'F#tx', ...
    'w-6dB', 'w-6dB', 'PSL', 'clutter', 'contrast', 'CNR', 'gCNR', 'SNR0', 'SNRgain', '10logN');
fprintf('%4s %6s %8s %6s %7s %7s %7s %7s %8s %9s %6s %6s %5s %8s %7s\n', '', '(deg)', '(mm)', '(dB)', '', ...
    'meas', 'theory', '(dB)', '(dB)', '(dB)', '', '', '', '(dB)', '(dB)');
% narrowband focal-plane prediction: receive pattern of the active aperture at z = 30 mm
% times the transmit factor |sum_n exp(j 2 pi x sin(theta_n)/lam)| / N
uz  = (xe - 0)/(30/(2*Fnum));  wz = double(abs(uz) <= 1);
if strcmp(rxApod, 'hann'), wz = wz .* (0.5 + 0.5*cos(pi*uz)); end
xn  = (0:0.002:3).';                                       % lateral offset (mm)
rxP = abs(exp(-1i*2*pi*xn*xe/(lam*30)) * wz.') / sum(wz);   % receive pattern, paraxial
for k = 1:nS
    idx = subsets{k};  N = numel(idx);
    sA = sum(imgA(:,idx), 2);  sB = sum(imgB(:,idx), 2);  sC = sum(imgC(:,idx), 2);
    eL = abs(sA(inL));  pk = max(eL);  prof(:,k) = 20*log10(eL/pk);
    eA = abs(sA(inB))/pk;  bA{k} = reshape(20*log10(eA + eps), size(XG));
    % -6 dB lateral width from the fine line (linear interpolation of the crossings)
    [~, i0] = max(eL);  y = eL/pk;
    iR = i0; while y(iR) > 0.5, iR = iR + 1; end
    iLf = i0; while y(iLf) > 0.5, iLf = iLf - 1; end
    xR = xl(iR-1) + (0.5 - y(iR-1))/(y(iR) - y(iR-1))*(xl(iR) - xl(iR-1));
    xLf = xl(iLf+1) - (0.5 - y(iLf+1))/(y(iLf) - y(iLf+1))*(xl(iLf+1) - xl(iLf));
    w6 = xR - xLf;
    PSL = 20*log10(max(eA(box & ~mlobe)));
    clut = 10*log10(sum(eA(box & ~mlobe).^2) / sum(eA(box & mlobe).^2));
    % cyst metrics on the envelope (contrast on intensity)
    eB = abs(sB(inB));  bB{k} = reshape(20*log10(eB/mean(eB(cout)) + eps), size(XG));
    ei = eB(cin); eo = eB(cout);
    contrast = 10*log10(mean(ei.^2)/mean(eo.^2));
    CNR  = abs(mean(ei) - mean(eo)) / sqrt(var(ei) + var(eo));
    edges = linspace(0, max([ei; eo]), 101);
    hi = histc(ei, edges)/numel(ei);  ho = histc(eo, edges)/numel(eo);
    gCNR = 1 - sum(min(hi, ho));
    SNR0 = mean(eo)/std(eo);
    snr  = pk^2 / mean(abs(sC(inB)).^2);
    if k == 1, snr1 = snr; end
    s = sin(theta(idx));  ds = mean(diff(s));
    if N > 1, step = mean(diff(theta(idx)))*180/pi; rep = lam/ds; Ftx = 1/(2*tan(max(abs(theta(idx)))));
    else,     step = NaN; rep = Inf; Ftx = Inf; end
    % ghost: highest level of the lateral profile around +-replica (if inside the line)
    ghost = NaN;
    if rep < 7.5
        gw = abs(abs(xl) - rep) <= 0.2 + 0.3*rep;  ghost = max(prof(gw, k));
    end
    txP = abs(exp(1i*2*pi*xn*s/lam) * ones(N,1)) / N;      % transmit factor
    P2  = rxP .* txP;  wTh = 2*xn(find(P2 < 0.5, 1));        % -6 dB full width (mm)
    res(k,:) = [N step rep ghost Ftx w6 wTh PSL clut contrast CNR gCNR SNR0 10*log10(snr/snr1)];
    fprintf('%4d %6.1f %8.2f %6.1f %7.2f %7.3f %7.3f %7.1f %8.1f %9.1f %6.2f %6.3f %5.2f %8.2f %7.2f\n', ...
        res(k,:), 10*log10(N));
end
fprintf('\ntheory, narrowband, rect: -6 dB width one-way 1.207 lam F# = %.3f mm, two-way 0.886 lam F# = %.3f mm\n', ...
    1.207*lam*Fnum, 0.886*lam*Fnum);
fprintf('F#tx = 1/(2 tan(theta_max)) = %.2f for +-%g deg; span needed for F#tx = F#rx = %.1f: +-%.1f deg\n', ...
    1/(2*tan(thetaMax_deg*pi/180)), thetaMax_deg, Fnum, atan(1/(2*Fnum))*180/pi);
fprintf('angle-step rule (replica outside a %.1f mm wide field): d(theta) <= lambda/W = %.2f deg\n', L, lam/L*180/pi);
fprintf('Montaldo count for focused-equivalent quality N = L/(lambda F#) = %.0f (L = %.1f mm, F# %.1f)\n', ...
    L/(lam*Fnum), L, Fnum);
fprintf('total runtime %.1f s\n', toc);

%% Figures
figure('Position', [50 50 1200 560]);
for k = 1:nS
    subplot(2, nS, k); imagesc(xg, zg, bA{k}, [-60 0]); axis image; colormap(gray);
    title(sprintf('%d angle(s): PSF, 60 dB', res(k,1))); xlabel('x (mm)'); ylabel('z (mm)');
    subplot(2, nS, nS+k); imagesc(xg, zg, bB{k}, [-40 10]); axis image; colormap(gray);
    title(sprintf('cyst: C %.1f dB, gCNR %.2f', res(k,10), res(k,12))); xlabel('x (mm)'); ylabel('z (mm)');
end
figure('Position', [80 80 1100 420]);
subplot(1,2,1); plot(xl, prof, 'LineWidth', 1.2); grid on; ylim([-70 5]);
xlabel('lateral position x (mm) at z = 30 mm'); ylabel('envelope (dB re peak)');
legend('1 angle', '3 angles (15 deg step)', '11 angles (3 deg step)', '31 angles (1 deg step)', 'Location', 'northeast');
title('Lateral profile: replicas at lambda/d(sin theta)');
subplot(1,2,2); Nv = res(:,1);
semilogx(Nv, res(:,14), 'o-', Nv, 10*log10(Nv), 'k--', Nv, res(:,8) - res(1,8), 's-', Nv, res(:,10) - res(1,10), 'd-', 'LineWidth', 1.2);
grid on; xlabel('number of angles N'); ylabel('dB (change re 1 angle)');
legend('SNR gain (measured)', '10 log_{10} N', 'peak sidelobe', 'cyst contrast', 'Location', 'southwest');
title('What compounding buys');
