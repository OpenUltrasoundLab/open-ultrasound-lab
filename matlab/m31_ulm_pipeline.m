%% m31_ulm_pipeline.m
% Purpose : a complete, minimal ULM pipeline on simulated beamformed IQ data:
%           synthesis -> SVD clutter filter -> matched-filter detection ->
%           sub-pixel localization (weighted average, Gaussian 3-point) ->
%           greedy frame-to-frame assignment -> track filter -> density map.
%           Every stage is compared with the ground truth that only a
%           simulation provides (true bubble positions and identities).
% Model   : 64 x 64 pixels of lambda/4 (5 MHz, lambda = 308 um), 1000 frames/s,
%           3 blocks of 100 frames. Six vessels (120, 60, 60, 30, 30, 20 um;
%           centreline speeds 25 ... 1.5 mm/s, parabolic profile). Bubbles:
%           Gaussian PSF, FWHM lambda (lateral) x 0.75 lambda (axial),
%           log-normal amplitude (4 dB), random phase, Doppler phase -4*pi*z/lambda.
%           Clutter: 3 speckle fields with slow complex gains, 30 dB above the
%           median bubble. White complex noise set by the median-bubble ENR.
%           Same model as Demo 31.1 (different random numbers).
% Change  : ENRdB, Nb (bubbles in the field), rList (SVD cut-offs), Pfa,
%           gate (linking distance, lambda), Lmin (minimum track length).
% Expect  : three clutter components (0-4 Hz mean Doppler, then a drop of
%           17 dB to bubble components at 50-100 Hz); r = 3 is the best cut-off
%           (recall 0.45, 0.26 false alarms per frame); r = 1 or 2 detects
%           almost nothing, r >= 4 loses bubbles and adds false alarms.
%           Localization error about 3x the median-bubble CRLB, mostly the
%           3-point estimator on lambda/4 pixels with white pixel noise (it
%           uses 3 samples of a PSF 4 pixels wide), also for isolated bubbles.
%           Greedy links 93 % correct. Across the 60 um vessel the localizations
%           spread 25.3 um while the true positions of the same (few) bubbles
%           spread 7.9 um: the quadrature difference, 24.0 um, matches the
%           actual rms perpendicular error (23.6 um). Runtime about 10 s (Octave).
rng(1);

%% Parameters
lam_um = 308; pix = 0.25; Nx = 64; Nz = 64; Np = Nx*Nz;  % pixel = lambda/4
FR = 1000; Tb = 100; nBlk = 3;                            % frames/s, frames per block
wx = 1/2.3548/pix; wz = 0.75/2.3548/pix;                  % PSF amplitude std (pixels)
ENRdB = 24; Nb = 20;                                      % median-bubble ENR, bubbles
Pfa = 1e-5; gate = 0.25; Lmin = 10;                       % detection, linking (lambda), track length
rList = [0 1 2 3 4 6 8 12];                               % SVD cut-offs to compare (block 1)
rUse = 3;                                                 % cut-off for the full pipeline
sigp = sqrt(pi*wx*wz/10^(ENRdB/10));                      % per-pixel complex noise std

%% Vessels: centrelines resampled to 0.02 lambda, normals, lengths
P = {[0.5 2;4 3;7 5.2;10 7.8;12.5 10.2;15.5 12.6], [7 5.2;5.6 7.8;4.6 10.2;3.6 12.6;3 15.5], ...
     [10 7.8;12 6;13.6 4.2;15.5 2.2], [4.6 10.2;6.4 11.8;8.2 13.6;9.4 15.5], ...
     [13.6 4.2;12 2.6;10.6 1.6;9 0.6], [12.5 10.2;11.6 12.2;12.2 14;12 15.5]};
dV = [120 60 60 30 30 20]; vV = [25 10 8 3 2.5 1.5];      % diameter (um), centreline speed (mm/s)
nV = numel(P); VX = cell(1,nV); VZ = VX; NXv = VX; NZv = VX; Lv = zeros(1,nV);
for k = 1:nV
    Q = P{k}; t = [0; cumsum(sqrt(sum(diff(Q).^2, 2)))];
    tt = linspace(0, t(end), 2000)';
    x = interp1(t, Q(:,1), tt, 'pchip'); z = interp1(t, Q(:,2), tt, 'pchip');
    cum = [0; cumsum(sqrt(diff(x).^2 + diff(z).^2))];
    su = (0:0.02:cum(end))';
    VX{k} = interp1(cum, x, su); VZ{k} = interp1(cum, z, su); Lv(k) = su(end);
    gx = gradient(VX{k}); gz = gradient(VZ{k}); gl = sqrt(gx.^2 + gz.^2);
    NXv{k} = -gz./gl; NZv{k} = gx./gl;
end
rV = dV/lam_um/2;                                         % radius (lambda)
vf = vV/(lam_um/1000)/FR;                                 % speed (lambda per frame)
wV = dV.*Lv; cw = cumsum(wV)/sum(wV);                     % bubbles per vessel ~ d x L (demo compromise)

%% Clutter fields (speckle = white noise convolved with the PSF)
kx = exp(-((-ceil(3*wx):ceil(3*wx)).^2)/(2*wx^2)); kz = exp(-((-ceil(3*wz):ceil(3*wz))'.^2)/(2*wz^2));
ampC = 10.^([30 18 6]/20); T = cell(1,3);
for k = 1:3
    W = conv2(kz, kx, randn(Nz,Nx) + 1i*randn(Nz,Nx), 'same');
    T{k} = W/sqrt(mean(abs(W(:)).^2))*ampC(k);
end
phc = 2*pi*rand(1,3); fC = [2.1 -5.3];                    % clutter gain phases, rotation (Hz)

%% Bubbles: vessel, arclength, lateral position u, amplitude, phase, id
bk = zeros(Nb,1); bs = bk; bu = bk; ba = bk; bph = bk; bid = (1:Nb)'; nid = Nb;
for i = 1:Nb
    bk(i) = find(rand <= cw, 1); bs(i) = rand*Lv(bk(i));
    u = 2*rand - 1; while rand > 1 - u^2, u = 2*rand - 1; end; bu(i) = u;
    ba(i) = 10^(4*randn/20); bph(i) = 2*pi*rand;
end

%% 1) Synthesis
tic;
nF = nBlk*Tb; IQ = zeros(Nz, Nx, nF); truth = cell(1, nF);
xg = 0:Nx-1; zg = (0:Nz-1)';
for f = 1:nF
    t = (f-1)/FR;
    c0 = exp(-1i*4*pi*0.02*sin(2*pi*3*t + phc(1)));
    fr = T{1}*c0 + T{2}*exp(1i*(2*pi*fC(1)*t + phc(2))) + T{3}*exp(1i*(2*pi*fC(2)*t + phc(3)));
    tl = zeros(0,4);
    for i = 1:Nb
        k = bk(i); bs(i) = bs(i) + vf(k)*(1 - bu(i)^2);
        if bs(i) > Lv(k)                                   % respawn at an inlet
            k = find(rand <= cw, 1); bk(i) = k; bs(i) = 0;
            u = 2*rand - 1; while rand > 1 - u^2, u = 2*rand - 1; end; bu(i) = u;
            ba(i) = 10^(4*randn/20); bph(i) = 2*pi*rand; nid = nid + 1; bid(i) = nid;
        end
        j = min(numel(VX{k}), max(1, round(bs(i)/0.02) + 1));
        x = VX{k}(j) + bu(i)*rV(k)*NXv{k}(j); z = VZ{k}(j) + bu(i)*rV(k)*NZv{k}(j);
        X = x/pix - 0.5; Z = z/pix - 0.5;                  % pixel coordinates (0-based)
        a = ba(i)*(1 + 0.1*randn)*exp(1i*(bph(i) - 4*pi*z));
        fr = fr + a*exp(-(zg - Z).^2/(2*wz^2))*exp(-(xg - X).^2/(2*wx^2));
        tl(end+1,:) = [x z bid(i) k];                      %#ok<SAGROW>
    end
    IQ(:,:,f) = fr + sigp/sqrt(2)*(randn(Nz,Nx) + 1i*randn(Nz,Nx));
    truth{f} = tl;
end
fprintf('Synthesis: %d frames of %d x %d pixels (%.1f s)\n', nF, Nz, Nx, toc);

%% 2) SVD of each block; singular spectrum and Doppler of the temporal vectors
U = cell(1,nBlk); Sg = U; Vv = U;
for b = 1:nBlk
    X = reshape(IQ(:,:,(b-1)*Tb + (1:Tb)), Np, Tb);
    [U{b}, Sg{b}, Vv{b}] = svd(X, 'econ');
end
sv = diag(Sg{1}); w = conj(Vv{1});                         % temporal signature of component k
fd = abs(angle(sum(w(2:end,:).*conj(w(1:end-1,:)), 1)))*FR/(2*pi);
fprintf('\nBlock 1: singular values (dB re sigma_1) and |mean Doppler| of components 1-8\n');
fprintf('  k      : %s\n', sprintf('%7d', 1:8));
fprintf('  dB     : %s\n', sprintf('%7.1f', 20*log10(sv(1:8)/sv(1))));
fprintf('  |fD| Hz: %s\n', sprintf('%7.1f', fd(1:8)));

%% 3-4) Detection and localization (anonymous helpers, base MATLAB)
kmx = exp(-((-ceil(2.5*wx):ceil(2.5*wx)).^2)/(2*wx^2)); kmz = exp(-((-ceil(2.5*wz):ceil(2.5*wz))'.^2)/(2*wz^2));
thr = log(1/Pfa);                                          % eta^2/sigma^2 for Pfa = exp(-eta^2/sigma^2)
filt = @(b, r) reshape(reshape(IQ(:,:,(b-1)*Tb + (1:Tb)), Np, Tb) - U{b}(:,1:r)*Sg{b}(1:r,1:r)*Vv{b}(:,1:r)', Nz, Nx, Tb);
o3 = @(a, b, c) max(-0.5, min(0.5, 0.5*(log(a) - log(c))/min(log(a) - 2*log(b) + log(c), -eps)));
res = zeros(numel(rList), 4);                              % recall, FA/frame, RMSE (um), median error (um)
for pass = 1:2
    if pass == 1, rr = rList; blocks = 1; else, rr = rUse; blocks = 1:nBlk; end
    for ir = 1:numel(rr)
        L = zeros(0, 9);                                   % [x z frame id ex ez xWA zWA vessel]
        nTrue = 0;
        for b = blocks
            F = filt(b, rr(ir));
            for t = 1:Tb
                f = (b-1)*Tb + t; A = abs(F(:,:,t));
                M = abs(conv2(kmz, kmx, F(:,:,t), 'same')).^2;
                s2 = median(M(:))/log(2);
                C = M(4:end-3, 4:end-3); pk = C > s2*thr;
                for di = -1:1, for dj = -1:1
                    if di || dj, pk = pk & C >= M((4:end-3)+di, (4:end-3)+dj); end
                end, end
                [ii, jj] = find(pk); ii = ii + 3; jj = jj + 3;
                tl = truth{f}; inner = tl(:,1) > 0.75 & tl(:,1) < 15.25 & tl(:,2) > 0.75 & tl(:,2) < 15.25;
                nTrue = nTrue + sum(inner);
                for q = 1:numel(ii)
                    blk = A(ii(q)-1:ii(q)+1, jj(q)-1:jj(q)+1); [~, m] = max(blk(:)); [di, dj] = ind2sub([3 3], m);
                    i0 = ii(q) + di - 2; j0 = jj(q) + dj - 2;
                    i0 = min(Nz-3, max(3, i0)); j0 = min(Nx-3, max(3, j0));
                    gxo = o3(A(i0,j0-1), A(i0,j0), A(i0,j0+1)); gzo = o3(A(i0-1,j0), A(i0,j0), A(i0+1,j0));
                    Wb = A(i0-2:i0+2, j0-2:j0+2); [jw, iw] = meshgrid(-2:2, -2:2);
                    xw = (j0 - 1 + sum(Wb(:).*jw(:))/sum(Wb(:)) + 0.5)*pix; zw = (i0 - 1 + sum(Wb(:).*iw(:))/sum(Wb(:)) + 0.5)*pix;
                    x = (j0 - 1 + gxo + 0.5)*pix; z = (i0 - 1 + gzo + 0.5)*pix;
                    d = sqrt((tl(:,1) - x).^2 + (tl(:,2) - z).^2); [dm, im] = min(d);
                    if dm < 0.5, L(end+1,:) = [x z f tl(im,3) x-tl(im,1) z-tl(im,2) xw zw tl(im,4)]; %#ok<SAGROW>
                    else, L(end+1,:) = [x z f -1 0 0 xw zw 0]; end %#ok<SAGROW>
                end
            end
        end
        tp = L(:,4) >= 0; [~, uq] = unique(L(tp,3)*1e6 + L(tp,4)); E = L(tp,:); E = E(uq,:);
        er = sqrt(E(:,5).^2 + E(:,6).^2)*lam_um;
        if pass == 1
            if isempty(er), er = NaN; end
            res(ir,:) = [size(E,1)/max(1,nTrue), sum(~tp)/(numel(blocks)*Tb), sqrt(mean(er.^2)), median(er)];
        end
    end
end
fprintf('\nSVD cut-off sweep on block 1 (Pfa = %.0e, Gaussian 3-point localizer)\n', Pfa);
fprintf('   r   recall   FA/frame   RMSE(um)   median(um)\n');
for ir = 1:numel(rList), fprintf('  %2d   %5.2f    %7.2f    %7.1f    %7.1f\n', rList(ir), res(ir,:)); end

crx = wx*pix/sqrt(10^(ENRdB/10))*lam_um; crz = wz*pix/sqrt(10^(ENRdB/10))*lam_um;
tp = L(:,4) >= 0;
eGx = L(tp,5)*lam_um; eGz = L(tp,6)*lam_um;
idx = find(tp); tl0 = zeros(numel(idx),2);
for q = 1:numel(idx), tl = truth{L(idx(q),3)}; m = find(tl(:,3) == L(idx(q),4), 1); tl0(q,:) = tl(m,1:2); end
eWx = (L(tp,7) - tl0(:,1))*lam_um; eWz = (L(tp,8) - tl0(:,2))*lam_um;
fprintf('\nFull pipeline, r = %d, %d blocks (%.1f s of data): %d localizations, %d true positives\n', rUse, nBlk, nF/FR, size(L,1), sum(tp));
fprintf('  CRLB of the median bubble (ENR %d dB): sigma_x = %.1f um, sigma_z = %.1f um\n', ENRdB, crx, crz);
fprintf('  Gaussian 3-point: RMSE x = %.1f um, z = %.1f um; median |error| = %.1f um\n', sqrt(mean(eGx.^2)), sqrt(mean(eGz.^2)), median(sqrt(eGx.^2 + eGz.^2)));
fprintf('  Weighted average: RMSE x = %.1f um, z = %.1f um; median |error| = %.1f um\n', sqrt(mean(eWx.^2)), sqrt(mean(eWz.^2)), median(sqrt(eWx.^2 + eWz.^2)));

%% Isolated versus crowded bubbles (nearest other true bubble beyond 1.5 lateral FWHM)
nn = zeros(numel(idx),1); amp = nn;
for q = 1:numel(idx)
    tl = truth{L(idx(q),3)}; m = tl(:,3) == L(idx(q),4);
    d = sqrt((tl(~m,1) - tl0(q,1)).^2 + ((tl(~m,2) - tl0(q,2))/0.75).^2); if isempty(d), d = Inf; end
    nn(q) = min(d);
end
iso = nn > 1.5; eG = sqrt(eGx.^2 + eGz.^2);
fprintf('  isolated (%d %%): RMSE x = %.1f, z = %.1f um, median |error| %.1f um; crowded: median |error| %.1f um\n', ...
    round(100*mean(iso)), sqrt(mean(eGx(iso).^2)), sqrt(mean(eGz(iso).^2)), median(eG(iso)), median(eG(~iso)));

%% Spread across the 60 um vessel (vessel 2) versus the prediction
% Lateral positions u are drawn with p(u) ~ 1 - u^2 (flux weighting of a 2D parabolic profile), so the
% model spread of the lumen is R/sqrt(5). The true positions of the matched bubbles give the actual
% lumen spread of this sample (a long steady-state run would drift toward uniform u, R/sqrt(3), because
% slow wall bubbles stay in the field longer; a cylindrical vessel gives R/sqrt(6) and R/2 instead).
in2 = tp & L(:,9) == 2; q2 = find(in2); perp = zeros(numel(q2),1); perpT = perp;
for q = 1:numel(q2)
    x = L(q2(q),1); z = L(q2(q),2); dd = (VX{2} - x).^2 + (VZ{2} - z).^2; [~, j] = min(dd);
    perp(q) = ((x - VX{2}(j))*NXv{2}(j) + (z - VZ{2}(j))*NZv{2}(j))*lam_um;
    m = find(idx == q2(q), 1); xt = tl0(m,1); zt = tl0(m,2);
    perpT(q) = ((xt - VX{2}(j))*NXv{2}(j) + (zt - VZ{2}(j))*NZv{2}(j))*lam_um;
end
sModel = dV(2)/2/sqrt(5);                                  % std of u*R with p(u) ~ 1 - u^2
fprintf('\nVessel 2 (60 um): %d localizations; spread across the vessel = %.1f um\n', numel(perp), std(perp));
fprintf('  true positions of the same bubbles: %.1f um (model R/sqrt(5) = %.1f um)\n', std(perpT), sModel);
fprintf('  implied perpendicular localization error = %.1f um (actual rms perpendicular error %.1f um)\n', ...
    sqrt(max(0, std(perp)^2 - std(perpT)^2)), sqrt(mean((perp - perpT).^2)));

%% 5) Tracking: greedy global assignment within the gate, then track filter
frames = L(:,3); ntr = 0; tid = zeros(size(L,1),1); nl = 0; nlt = 0; nlc = 0;
for f = 1:nF
    cur = find(frames == f); if isempty(cur), continue; end
    prv = find(frames == f-1);
    if ~isempty(prv)
        D = sqrt((L(prv,1) - L(cur,1)').^2 + (L(prv,2) - L(cur,2)').^2);
        [dsort, ord] = sort(D(:)); usedA = false(numel(prv),1); usedB = false(numel(cur),1);
        for q = 1:numel(ord)
            if dsort(q) > gate, break; end
            [a, b2] = ind2sub(size(D), ord(q));
            if ~usedA(a) && ~usedB(b2)
                usedA(a) = true; usedB(b2) = true; tid(cur(b2)) = tid(prv(a)); nl = nl + 1;
                if L(prv(a),4) >= 0 && L(cur(b2),4) >= 0, nlt = nlt + 1; nlc = nlc + (L(prv(a),4) == L(cur(b2),4)); end
            end
        end
    end
    for q = cur', if tid(q) == 0, ntr = ntr + 1; tid(q) = ntr; end, end
end
len = accumarray(tid, 1); keep = find(len >= Lmin);
fprintf('\nTracking (greedy, gate %.2f lambda = %.0f um): %d links, %.1f %% of true-true links correct\n', gate, gate*lam_um, nl, 100*nlc/max(1,nlt));
fprintf('  %d tracks, %d kept (>= %d frames), %.0f %% of localizations kept\n', ntr, numel(keep), Lmin, 100*sum(len(keep))/size(L,1));

%% 6) Render (lambda/12 grid), each track counted once per pixel
RG = 192; rp = 16/RG; den = zeros(RG);
for k = keep'
    q = find(tid == k); [~, o] = sort(L(q,3)); q = q(o); st = false(RG);
    for m = 2:numel(q)
        n = max(1, ceil(hypot(L(q(m),1) - L(q(m-1),1), L(q(m),2) - L(q(m-1),2))/(0.3*rp)));
        xs = L(q(m-1),1) + (0:n)/n*(L(q(m),1) - L(q(m-1),1)); zs = L(q(m-1),2) + (0:n)/n*(L(q(m),2) - L(q(m-1),2));
        ix = floor(xs/rp) + 1; iz = floor(zs/rp) + 1; ok = ix >= 1 & ix <= RG & iz >= 1 & iz <= RG;
        st(sub2ind([RG RG], iz(ok), ix(ok))) = true;
    end
    den = den + st;
end
mask = false(RG); [gxr, gzr] = meshgrid(((1:RG) - 0.5)*rp, ((1:RG) - 0.5)*rp);
for k = 1:nV
    for j = 1:4:numel(VX{k})
        mask = mask | ((gxr - VX{k}(j)).^2 + (gzr - VZ{k}(j)).^2 <= (max(rV(k), 0.5*rp) + 0.5*rp)^2);
    end
end
fprintf('  saturation: %.1f %% of true-lumen pixels visited by at least one kept track\n', 100*mean(den(mask) > 0));
fprintf('Total runtime %.1f s\n', toc);

%% Figures
figure('Position', [50 50 1100 650]);
subplot(2,3,1); imagesc((0:Nx-1)*pix*lam_um/1000, (0:Nz-1)*pix*lam_um/1000, 20*log10(abs(IQ(:,:,50))/max(max(abs(IQ(:,:,50)))))); caxis([-40 0]);
colormap(gray); axis image; title('raw IQ frame 50 (dB)'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,2); stem(1:20, 20*log10(sv(1:20)/sv(1)), 'filled'); hold on; plot([rUse rUse] + 0.5, [-70 0], 'r--');
xlabel('component k'); ylabel('\sigma_k (dB re \sigma_1)'); title('block 1 singular values'); ylim([-70 5]);
Ff = filt(1, rUse); subplot(2,3,3); imagesc((0:Nx-1)*pix*lam_um/1000, (0:Nz-1)*pix*lam_um/1000, 20*log10(abs(Ff(:,:,50))/max(max(abs(Ff(:,:,50)))))); caxis([-30 0]);
axis image; hold on; tl = truth{50}; plot(tl(:,1)*lam_um/1000, tl(:,2)*lam_um/1000, 'r+'); q = L(:,3) == 50; plot(L(q,1)*lam_um/1000, L(q,2)*lam_um/1000, 'co');
title('filtered frame 50, detections (o), truth (+)'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,4); plot(L(:,1)*lam_um/1000, L(:,2)*lam_um/1000, 'b.', 'MarkerSize', 2); axis ij; axis image; axis([0 4.93 0 4.93]);
title('all localizations'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,5); hold on; for k = keep', q = find(tid == k); plot(L(q,1)*lam_um/1000, L(q,2)*lam_um/1000, '-'); end
axis ij; axis image; axis([0 4.93 0 4.93]); box on; title(sprintf('%d kept tracks', numel(keep))); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,6); imagesc(((1:RG) - 0.5)*rp*lam_um/1000, ((1:RG) - 0.5)*rp*lam_um/1000, log(1 + den)); axis image; colormap(gca, hot);
title('density map (log), \lambda/12 grid'); xlabel('x (mm)'); ylabel('z (mm)');
