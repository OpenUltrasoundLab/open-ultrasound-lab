%% m06_svd_clutter_filter.m  --  Module 06: SVD clutter filter on a simulated ultrafast ensemble
% Purpose : Simulate an ultrafast IQ ensemble (pixels x frames): moving tissue
%           speckle, a slow-flow vessel, a faster vessel, a few microbubble tracks,
%           white noise. Build the Casorati matrix, take its SVD, look at the
%           singular-value curve (with the Marchenko-Pastur noise edge), the
%           Doppler spectra of the temporal singular vectors, and power Doppler
%           for: raw data, an SVD band, and three FIXED temporal filters (first
%           difference, 2nd-order Butterworth IIR high-pass run forward-backward,
%           polynomial-regression filter). Then sweep tissue velocity and watch
%           the best SVD band move while every fixed temporal filter collapses.
%           Last, check the fast "Gram route" (eig of the n x n matrix S'S) in
%           double and single precision.
% Change  : vt_main (tissue peak axial velocity), the sweep list vts, PRF, n (frames),
%           blood_dB, noise_dB, the cutoffs fcs and polynomial degrees pdeg.
% Expect  : (rng(1)) at vt = 0 every filter works (SVD 24.7 dB, IIR 23.1 dB contrast).
%           At vt = 8 mm/s (52 Hz peak tissue Doppler) the SVD band 9-46 gives
%           15.2 dB vessel/background contrast; IIR, regression and first difference
%           all stay near 0 dB. The oracle k_lo rises 2 -> 5 -> 9 -> 11 -> 15 for
%           vt = 0, 4, 8, 16, 24 mm/s; k_hi from the MP edge is 45-47. The Gram
%           route matches svd() to ~1e-11 in double; in single it misplaces
%           sigma_k by > 1 % about 66-69 dB below sigma_1.
%           These are oracle (ground-truth-mask) choices: in vivo you have no masks.
% Runtime : ~20 s in Octave (five simulated ensembles; the simulation dominates).
clear; close all; rng(1);

%% 1. Parameters
c0 = 1540; f0 = 5e6; lam = c0/f0*1e3;        % wavelength (mm), 0.308 mm
k2 = 4*pi/lam;                               % two-way wavenumber (rad/mm)
PRF = 1000; n = 128;                         % frame rate (Hz), frames -> 128 ms
NX = 48; NZ = 48; dx = 0.1;                  % pixels, pixel pitch (mm)
xg = (0:NX-1)*dx; zg = (0:NZ-1)*dx; m = NX*NZ;
sx = 0.17; sz = 0.10;                        % PSF std lateral / axial (mm)
blood_dB = -30;                              % blood scatterer amplitude re tissue scatterer
noise_dB = -50;                              % noise power re mean tissue pixel power
fm = 3;                                      % tissue motion frequency (Hz), e.g. breathing/pulsation
vt_main = 8;                                 % tissue peak axial velocity for the main figure (mm/s)
vts = [0 4 8 16 24];                         % velocity sweep (mm/s)
fcs = [25 50 100 150];                       % high-pass cutoffs to try (Hz)
pdeg = [1 2 4 8];                            % polynomial-regression degrees to try

% Scatterers are drawn once and reused for every velocity, so only motion changes.
L = NX*dx; Ntis = 4000;
tis = [rand(Ntis,1)*(L+2)-1, rand(Ntis,1)*(L+2)-1, (randn(Ntis,1)+1i*randn(Ntis,1))/sqrt(2)];
ves = [1.2 2.2 10 0.30  5;                    % [x0 z0 angle(deg) radius(mm) vmax(mm/s)]: slow vessel
       3.4 2.6 70 0.22 15];                   % faster, steep vessel
blood = [];                                   % [vessel, s0, r, amplitude]
for v = 1:size(ves,1)
    nb = round(500*ves(v,4)/0.3);
    blood = [blood; v*ones(nb,1), rand(nb,1)*8-4, (2*rand(nb,1)-1)*ves(v,4), ...
             10^(blood_dB/20)*(randn(nb,1)+1i*randn(nb,1))/sqrt(2)]; %#ok<AGROW>
end
% microbubbles: [xc zc radius(mm) speed(mm/s) start angle(rad)] on circular arcs
mb = [2.4 1.2 1.5 12 0.3; 1.6 3.6 1.2 20 2.0; 3.8 3.8 1.8 9 4.0; 0.8 1.0 1.0 16 5.5];
mbAmp = 1.0;                                  % 0 dB re a tissue scatterer

PSFx = @(x) exp(-(xg(:) - x(:).').^2/(2*sx^2));        % NX x Ns
PSFz = @(z) exp(-(zg(:) - z(:).').^2/(2*sz^2));        % NZ x Ns

%% 2. Simulate one ensemble per tissue velocity; each gives a Casorati matrix S (m x n)
S_all = cell(numel(vts),1); masks = [];
for iv = 1:numel(vts)
    vt = vts(iv); w = 2*pi*fm;
    S = zeros(m, n);
    for f = 1:n
        t = (f-1)/PRF;
        dz0 = vt/w*sin(w*t);                  % bulk axial displacement (mm), peak velocity vt
        % mild non-rigid motion: displacement grows 30 % across the field (depth-dependent)
        dzs = @(z) dz0*(1 + 0.3*(z - L/2)/L);
        dxs = 0.35*dz0;
        zt = tis(:,2) + dzs(tis(:,2)); xt = tis(:,1) + dxs;
        amp = tis(:,3).*exp(-1i*k2*zt);
        img = PSFx(xt)*(amp.*PSFz(zt).');       % NX x NZ
        % blood: parabolic profile, riding on the tissue motion
        a = ves(blood(:,1),3)*pi/180; ux = cos(a); uz = sin(a); R = ves(blood(:,1),4);
        vel = ves(blood(:,1),5).*(1 - (blood(:,3)./R).^2);
        sp = mod(blood(:,2) + vel*t + 4, 8) - 4;
        xb = ves(blood(:,1),1) + ux.*sp - uz.*blood(:,3);
        zb = ves(blood(:,1),2) + uz.*sp + ux.*blood(:,3);
        zb = zb + dzs(zb); xb = xb + dxs;
        img = img + PSFx(xb)*((blood(:,4).*exp(-1i*k2*zb)).*PSFz(zb).');
        % microbubbles on arcs
        ang = mb(:,5) + mb(:,4)./mb(:,3)*t;
        xm = mb(:,1) + mb(:,3).*cos(ang) + dxs; zm = mb(:,2) + mb(:,3).*sin(ang);
        zm = zm + dzs(zm);
        img = img + PSFx(xm)*((mbAmp*exp(-1i*k2*zm)).*PSFz(zm).');
        S(:, f) = img(:);                     % column = one frame (x-major, NX*NZ)
    end
    if iv == 1, p0 = mean(abs(S(:,1)).^2); end % mean tissue pixel power (static case)
    sn = sqrt(p0*10^(noise_dB/10));           % noise std per entry
    S = S + sn*(randn(m,n) + 1i*randn(m,n))/sqrt(2);
    S_all{iv} = S;
end

%% 3. Masks (ground truth, at the ensemble-mean position; only used for scoring)
[XX, ZZ] = ndgrid(xg, zg);                    % NX x NZ, matches img(:)
dn = inf(NX, NZ); dv = cell(1, size(ves,1));
for v = 1:size(ves,1)
    a = ves(v,3)*pi/180;
    dv{v} = abs(-(XX - ves(v,1))*sin(a) + (ZZ - ves(v,2))*cos(a))/ves(v,4);   % distance in radii
    dn = min(dn, dv{v});
end
dmb = inf(NX, NZ); tt = linspace(0, (n-1)/PRF, 60);
for b = 1:size(mb,1)
    ang = mb(b,5) + mb(b,4)/mb(b,3)*tt;
    for q = 1:numel(tt)
        dmb = min(dmb, hypot(XX - mb(b,1) - mb(b,3)*cos(ang(q)), ZZ - mb(b,2) - mb(b,3)*sin(ang(q))));
    end
end
vesM = dn(:) < 0.6;                           % both vessel cores
slowM = dv{1}(:) < 0.6; fastM = dv{2}(:) < 0.6;
bgM = dn(:) > 2.5 & dmb(:) > 0.5;             % background, away from vessels and bubble tracks
contrast = @(pd) 10*log10(mean(pd(vesM))/mean(pd(bgM)));
cSlow = @(pd) 10*log10(mean(pd(slowM))/mean(pd(bgM)));
cFast = @(pd) 10*log10(mean(pd(fastM))/mean(pd(bgM)));
c3 = @(pd) [contrast(pd), cSlow(pd), cFast(pd)];

%% 4. Filters
pdOf = @(Y) mean(abs(Y).^2, 2);               % power Doppler per pixel
fd = @(Y) Y(:,2:end) - Y(:,1:end-1);          % first difference: |H| = 2|sin(pi f/PRF)|
% 2nd-order Butterworth high-pass via the bilinear transform (no toolbox)
bwb = @(fc) [1, -2, 1] / (1 + sqrt(2)*tan(pi*fc/PRF) + tan(pi*fc/PRF)^2);
bwa = @(fc) [1, 2*(tan(pi*fc/PRF)^2 - 1), 1 - sqrt(2)*tan(pi*fc/PRF) + tan(pi*fc/PRF)^2] ...
            / (1 + sqrt(2)*tan(pi*fc/PRF) + tan(pi*fc/PRF)^2);
hp1 = @(Y, fc) filter(bwb(fc), bwa(fc), Y - Y(:,1), [], 2);   % one pass; subtract first sample to limit the step
rev = @(Y) Y(:, end:-1:1);
iir = @(Y, fc) rev(hp1(rev(hp1(Y, fc)), fc));                 % forward-backward: zero phase, |H|^2 (like filtfilt)
cols = @(Y, idx) Y(:, idx);
iirPD = @(Y, fc) mean(abs(cols(iir(Y, fc), 21:n-20)).^2, 2);  % drop 20-frame transients at both ends
% Polynomial regression filter: remove each pixel's projection onto polynomials of
% degree 0..p in slow time. Like the SVD filter it is a projection, S - S*Q*Q',
% but onto a FIXED temporal subspace chosen before seeing the data.
tn = linspace(-1, 1, n)';                     % normalized slow time
regPD = @(Y, p) pdOf(Y - (Y*orth(tn.^(0:p)))*orth(tn.^(0:p))');

%% 5. Main case: singular values, Doppler of temporal singular vectors, PD images
iv = find(vts == vt_main); S = S_all{iv};
[U, Sg, V] = svd(S, 'econ'); s = diag(Sg);
sn = sqrt(p0*10^(noise_dB/10));
edgeMP = sn*(sqrt(m) + sqrt(n));              % Marchenko-Pastur upper edge of the noise bulk
k_hi = find(s > 1.05*edgeMP, 1, 'last');      % keep everything above the noise bulk
fprintf('Casorati %d x %d, vt = %g mm/s, tissue peak Doppler 2vt/lambda = %.1f Hz\n', m, n, vt_main, 2*vt_main/lam);
fprintf('sigma_k (dB re sigma_1), k = 1..10: %s\n', sprintf('%.0f ', 20*log10(s(1:10)/s(1))));
fprintf('MP noise edge = %.1f dB re sigma_1  ->  k_hi = %d (measured last sigma = %.1f dB)\n', ...
    20*log10(edgeMP/s(1)), k_hi, 20*log10(s(end)/s(1)));
NF = 256; fax = ((0:NF-1) - NF/2)*PRF/NF;
Pk = abs(fftshift(fft(conj(V), NF, 1), 1)).^2;           % NF x n, spectrum of each temporal vector
fbar = (abs(fax)*Pk)./sum(Pk, 1);                         % mean |Doppler| per rank (Hz)
fprintf('mean |Doppler| of temporal singular vectors, k = 1..12 (Hz): %s\n', sprintf('%.0f ', fbar(1:12)));
% PD for every band start k_lo (with k_hi fixed): uses PD = sum sigma^2 |u|^2 / n
E = (abs(U).^2) .* (s.'.^2) / n;              % m x n, energy of each rank per pixel
conK = zeros(1, k_hi);
for klo = 1:k_hi, conK(klo) = contrast(sum(E(:, klo:k_hi), 2)); end
[cbest, klo_best] = max(conK);
pdRaw = pdOf(S); pdSVD = sum(E(:, klo_best:k_hi), 2);
pdFD = pdOf(fd(S));
cI3 = zeros(numel(fcs), 3); pdI = cell(1, numel(fcs));
for j = 1:numel(fcs), pdI{j} = iirPD(S, fcs(j)); cI3(j,:) = c3(pdI{j}); end
[cIb, jI] = max(cI3(:,1)); pdIIR = pdI{jI};
fprintf('\nPower Doppler contrast re background (dB):   both   slow   fast vessel\n');
fprintf('  raw (no filter)                        %6.1f %6.1f %6.1f\n', c3(pdRaw));
fprintf('  SVD band %2d-%2d (oracle k_lo)           %6.1f %6.1f %6.1f\n', klo_best, k_hi, c3(pdSVD));
fprintf('  first difference                       %6.1f %6.1f %6.1f\n', c3(pdFD));
for j = 1:numel(fcs)
    fprintf('  IIR Butterworth HP, zero-phase %3d Hz  %6.1f %6.1f %6.1f\n', fcs(j), cI3(j,:));
end
for j = 1:numel(pdeg)
    fprintf('  polynomial regression, degree %d        %6.1f %6.1f %6.1f\n', pdeg(j), c3(regPD(S, pdeg(j))));
end
fprintf('Doppler resolution of the ensemble PRF/n = %.1f Hz (%.2f mm/s axial)\n', PRF/n, PRF/n*lam/2);

%% 6. Velocity sweep: the best band moves, the high-pass fails
fprintf('\nContrast both / slow / fast vessel (dB) versus tissue velocity\n');
fprintf(' vt(mm/s) 2vt/lam(Hz)  oracle band   SVD band             best IIR (fc)              best regr. (p)   1st diff\n');
conCurves = zeros(numel(vts), 40);
for iv = 1:numel(vts)
    Sv = S_all{iv}; [Uv, Sv2, ~] = svd(Sv, 'econ'); sv = diag(Sv2);
    kh = find(sv > 1.05*edgeMP, 1, 'last');
    Ev = (abs(Uv).^2) .* (sv.'.^2) / n;
    cc = -inf(1, kh); for klo = 1:kh, cc(klo) = contrast(sum(Ev(:, klo:kh), 2)); end
    [cb, kb] = max(cc); conCurves(iv, 1:min(40,kh)) = cc(1:min(40,kh));
    ci = -inf; for j = 1:numel(fcs), pj = iirPD(Sv, fcs(j)); cj = contrast(pj); if cj > ci, ci = cj; jb = j; pib = pj; end, end
    cr = -inf; for j = 1:numel(pdeg), cj = contrast(regPD(Sv, pdeg(j))); if cj > cr, cr = cj; jr = j; end, end
    fprintf('  %5.0f     %6.1f       %2d-%2d   %5.1f %5.1f %5.1f    %5.1f %5.1f %5.1f (%3d Hz)    %5.1f (%d)       %5.1f\n', vts(iv), 2*vts(iv)/lam, kb, kh, ...
        c3(sum(Ev(:, kb:kh), 2)), c3(pib), fcs(jb), cr, pdeg(jr), contrast(pdOf(fd(Sv))));
end

%% 7. The fast route: eigendecomposition of the n x n Gram matrix S'S
% sigma_k^2 = eig(S'S), V = its eigenvectors; the filter S*Vb*Vb' never needs U.
% Forming S'S squares the condition number: a sigma_k below sqrt(eps)*sigma_1 is lost.
S = S_all{find(vts == vt_main)};
G = S'*S; G = (G + G')/2;                     % n x n Hermitian, cost O(m n^2)
[W, L] = eig(G); [lg, o] = sort(real(diag(L)), 'descend'); W = W(:, o);
sgD = sqrt(max(lg, 0));
Gs = single(S)'*single(S); Gs = (Gs + Gs')/2;
sgS = sqrt(max(double(sort(real(eig(Gs)), 'descend')), 0));
reS = abs(sgS - s)./s; kbad = find(reS > 0.01, 1);
Vb = W(:, klo_best:k_hi); pdGram = sum(abs(S*Vb).^2, 2)/n;   % = (1/n) sum_t |S Vb Vb'|^2, since Vb'Vb = I
fprintf('\nGram route (eig of %dx%d S''S): max rel. error of sigma_1..k_hi = %.1e (double)\n', n, n, max(abs(sgD(1:k_hi) - s(1:k_hi))./s(1:k_hi)));
fprintf('  PD from Gram route vs svd: max rel. difference %.1e\n', max(abs(pdGram - pdSVD)./pdSVD));
fprintf('  single precision: first sigma with > 1 %% error is k = %d at %.1f dB re sigma_1 (sqrt(eps_single) = %.1f dB)\n', ...
    kbad, 20*log10(s(kbad)/s(1)), 20*log10(sqrt(eps('single'))));

%% 8. Figures
figure('Position', [40 40 1100 640]);
subplot(2,3,1);
plot(20*log10(s/s(1)), 'k.-'); hold on; grid on;
plot([1 n], 20*log10(edgeMP/s(1))*[1 1], 'b--');
plot([klo_best klo_best], [-90 0], 'r'); plot([k_hi k_hi], [-90 0], 'r');
xlabel('singular value index k'); ylabel('\sigma_k (dB re \sigma_1)'); ylim([-90 5]);
title(sprintf('Singular values, vt = %g mm/s', vt_main)); legend('\sigma_k', 'MP noise edge', 'band', 'Location', 'northeast');
subplot(2,3,2);
Pn = 10*log10(Pk./max(Pk, [], 1));
imagesc(1:n, fax, Pn, [-30 0]); axis xy; colormap(gray); hold on;
plot([1 n], 2*vt_main/lam*[1 1], 'r--'); plot([1 n], -2*vt_main/lam*[1 1], 'r--');
xlabel('rank k'); ylabel('Doppler frequency (Hz)'); title('Spectra of temporal singular vectors (dB)');
subplot(2,3,3);
for iv = 1:numel(vts)
    kk = find(isfinite(conCurves(iv,:)) & conCurves(iv,:) ~= 0);
    plot(kk, conCurves(iv,kk), 'LineWidth', 1.5); hold on;
end
grid on; xlabel('k_{lo} (k_{hi} from MP edge)'); ylabel('contrast (dB)'); xlim([1 30]);
title('Best k_{lo} moves with tissue velocity');
legend(arrayfun(@(v) sprintf('%g mm/s', v), vts, 'UniformOutput', false), 'Location', 'southeast');
imgs = {pdRaw, pdSVD, pdIIR};
ttl = {'raw PD', sprintf('SVD band %d-%d', klo_best, k_hi), sprintf('IIR high-pass %d Hz', fcs(jI))};
for j = 1:3
    subplot(2,3,3+j);
    P = reshape(imgs{j}, NX, NZ).';           % rows = depth
    imagesc(xg, zg, 10*log10(P/max(P(:))), [-35 0]); axis image;
    xlabel('lateral (mm)'); ylabel('depth (mm)');
    title(sprintf('%s: %.1f dB', ttl{j}, contrast(imgs{j})));
end
