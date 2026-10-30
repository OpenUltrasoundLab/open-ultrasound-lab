%% m34_display_chain.m  --  from envelope to luminance: compression, gray map, bits, monitor, ambient
% Purpose : One synthetic B-mode frame (3 MHz, 60 % bandwidth, F#2, 40 x 40 mm) with an anechoic
%           cyst, a -6 dB lesion, three point targets, 0.5 dB/cm/MHz attenuation and a thermal
%           noise floor is pushed through the display chain:
%             envelope A -> TGC -> auto gain (median tissue at gray fraction g_t)
%             -> log compression with dynamic range DR, or mu-law (mu = 255)
%             -> gray map (linear, S-curve, reject) -> b-bit P-values
%             -> monitor: GSDF-calibrated or gamma 2.2 between Lmin and Lmax, plus ambient L_amb
%             -> JND index j(L') (DICOM PS3.14 inverse polynomial).
%           Metrics at every stage: CNR = |mu_in - mu_out| / sqrt(s_in^2 + s_out^2) and
%           gCNR = 1 - overlap (quantile bins of the pooled data, as in m26_cyst_cnr_gcnr.m),
%           for the anechoic cyst and the -6 dB lesion; distinct P-values used; display levels
%           separated by at least 1 JND; lesion mean difference in JNDs.
% Change  : cases table (DR, compression, map, bits, monitor, ambient), att, noise_dB, g_t.
% Expect  : gCNR identical on the envelope and on any strictly monotone, unclipped transform;
%           it drops only where clipping (DR, reject) or coarse quantization merges values.
%           CNR changes with every map. A gamma-2.2 monitor or bright ambient light shrinks the
%           lesion's mean difference in JNDs, not its gCNR in P-values.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(7);
c = 1540;  f0 = 3e6;  B = 0.6*f0;  Fn = 2;  lam = c/f0;
Nz = 256;  Nx = 256;  dz = 40e-3/Nz;  dx = 40e-3/Nx;  z0 = 10e-3;   % image starts at 10 mm depth
zc = z0 + ((0:Nz-1)' + 0.5)*dz;  xc = ((0:Nx-1) - Nx/2 + 0.5)*dx;
[X, Z] = meshgrid(xc, zc);

%% pulse-echo transfer function (separable, focal-plane approximation everywhere)
uz = [0:Nz/2-1, -Nz/2:-1]'/(Nz*dz);  ux = [0:Nx/2-1, -Nx/2:-1]/(Nx*dx);
sf = B/(2*sqrt(2*log(2)));
H = exp(-(c*uz/2).^2/(2*sf^2)) * max(0, 1 - abs(ux)*lam*Fn);

%% reflectivity: background 1, anechoic cyst, -6 dB lesion, point targets
Rc = 4e-3;
r1 = hypot(X - (-8e-3), Z - 25e-3);   r2 = hypot(X - 8e-3, Z - 32e-3);
refl = ones(Nz, Nx);  refl(r1 < Rc) = 0;  refl(r2 < Rc) = 10^(-6/10);
W = (randn(Nz, Nx) + 1i*randn(Nz, Nx))/sqrt(2) .* sqrt(refl);
h = ifft2(H);  apt = sqrt(1000*sum(abs(h(:)).^2)/max(abs(h(:)).^2));   % point peak 30 dB over speckle
pts = [15e-3 0; 25e-3 6e-3; 42e-3 -4e-3];
for q = 1:size(pts, 1)
    [~, iz] = min(abs(zc - pts(q, 1)));  [~, ix] = min(abs(xc - pts(q, 2)));
    W(iz, ix) = W(iz, ix) + apt;
end
zf = ifft2(fft2(W) .* H);
zf = zf / sqrt(mean(abs(zf(refl == 1)).^2));                 % unit mean speckle intensity
att = 0.5;  noise_dB = -48;                                  % dB/cm/MHz; noise re speckle at the surface
attdB = 2*att*(f0/1e6)*(zc*100);                             % two-way attenuation (dB)
zr = bsxfun(@times, zf, 10.^(-attdB/20)) + sqrt(10^(noise_dB/10)/2)*(randn(Nz, Nx) + 1i*randn(Nz, Nx));
A = bsxfun(@times, abs(zr), 10.^(attdB/20));                 % envelope after full TGC
fprintf('two-way attenuation %.0f to %.0f dB over the image; after TGC the noise floor rises to %.1f dB re speckle at the bottom\n', ...
        attdB(1), attdB(end), noise_dB + attdB(end));

%% ROIs
inC = r1 < 0.8*Rc;   outC = r1 > 1.3*Rc & r1 < 1.9*Rc;
inL = r2 < 0.8*Rc;   outL = r2 > 1.3*Rc & r2 < 1.9*Rc;
tissue = refl == 1 & ~(r1 < 1.3*Rc) & ~(r2 < 1.3*Rc);

%% helpers (anonymous only)
cnr   = @(v, a, b) abs(mean(v(a)) - mean(v(b))) / sqrt(var(v(a)) + var(v(b)));
qedge = @(p, nb) unique(p(max(1, round(linspace(1, numel(p), nb + 1)))));
fold  = @(hh) [hh(1:end-2); hh(end-1) + hh(end)];
ovl   = @(x, y, ed) sum(min(fold(histc(x(:), ed))/numel(x), fold(histc(y(:), ed))/numel(y)));
gcnrq = @(x, y) 1 - ovl(x, y, qedge(sort([x(:); y(:)]), 64));
gc    = @(v, a, b) gcnrq(v(a), v(b));
a = -1.3011877;   b = -2.5840191e-2; c1 = 8.0242636e-2; d = -1.0320229e-1; e = 1.3646699e-1;
f = 2.8745620e-2; g = -2.5468404e-2; hh = -3.1978977e-3; k = 1.2992634e-4; m = 1.3635334e-3;
Lj = @(j) 10.^((a + c1*log(j) + e*log(j).^2 + g*log(j).^3 + m*log(j).^4) ./ ...
             (1 + b*log(j) + d*log(j).^2 + f*log(j).^3 + hh*log(j).^4 + k*log(j).^5));
jL = @(L) 71.498068 + 94.593053*log10(L) + 41.912053*log10(L).^2 + 9.8247004*log10(L).^3 + ...
          0.28175407*log10(L).^4 - 1.1878455*log10(L).^5 - 0.18014349*log10(L).^6 + ...
          0.14710899*log10(L).^7 - 0.017046845*log10(L).^8;
Smap = @(x) (tanh(4*(x - 0.5)) - tanh(-2)) / (tanh(2) - tanh(-2));   % S-curve, 0 -> 0, 1 -> 1

%% chain: one row per case
% {compression, DR (dB), gray map, bits, monitor, ambient (cd/m^2), ambient assumed at calibration}
cases = { 'log', 60, 'linear', 8,  'gsdf',  0.2, 0.2;
          'log', 40, 'linear', 8,  'gsdf',  0.2, 0.2;
          'log', 80, 'linear', 8,  'gsdf',  0.2, 0.2;
          'log', 20, 'linear', 8,  'gsdf',  0.2, 0.2;
          'mu',  48, 'linear', 8,  'gsdf',  0.2, 0.2;
          'log', 60, 'S',      8,  'gsdf',  0.2, 0.2;
          'log', 60, 'reject', 8,  'gsdf',  0.2, 0.2;
          'log', 60, 'linear', 4,  'gsdf',  0.2, 0.2;
          'log', 60, 'linear', 10, 'gsdf',  0.2, 0.2;
          'log', 60, 'linear', 8,  'gamma', 0.2, 0.2;
          'log', 60, 'linear', 8,  'gsdf',  2,   2;
          'log', 60, 'linear', 8,  'gsdf',  50,  50;
          'log', 60, 'linear', 8,  'gsdf',  50,  0.2;
          'log', 60, 'linear', 8,  'gamma', 50,  0.2 };
Lmin = 0.5;  Lmax = 400;  g_t = 0.55;
D50 = median(20*log10(A(tissue)));
fprintf('\nenvelope (linear): CNR cyst %.3f lesion %.3f | gCNR cyst %.3f lesion %.3f\n', ...
        cnr(A, inC, outC), cnr(A, inL, outL), gc(A, inC, outC), gc(A, inL, outL));
fprintf('\n comp  DR  map     b  monitor amb (cal) | CNR cyst/lesion: gray    JND   | gCNR cyst/lesion: gray     JND    | P used  >=1JND  JNDs spanned  lesion dJ\n');
for q = 1:size(cases, 1)
    [cmp, DR, mp, nb, mon, amb, ambc] = deal(cases{q, :});
    Aref = 10^((D50 + DR*(1 - g_t))/20);                     % auto gain: median tissue at gray g_t
    if strcmp(cmp, 'log')
        x = min(1, max(0, (20*log10(A/Aref) + DR)/DR));
    else
        mu = 255;  x = log(1 + mu*min(1, A/Aref)) / log(1 + mu);
    end
    switch mp
        case 'linear', y = x;
        case 'S',      y = Smap(x);
        case 'reject', y = max(0, (x - 0.15)/0.85);
    end
    pmax = 2^nb - 1;  P = round(pmax*y);
    if strcmp(mon, 'gsdf')                                   % calibrated with ambient ambc included
        j0 = jL(Lmin + ambc);  j1 = jL(Lmax + ambc);
        Lpanel = Lj(j0 + P/pmax*(j1 - j0)) - ambc;
    else
        Lpanel = Lmin + (Lmax - Lmin)*(P/pmax).^2.2;
    end
    Lv = Lpanel + amb;  J = jL(Lv);
    lev = unique(J(:));  nd = 1;  last = lev(1);
    for t = 2:numel(lev), if lev(t) - last >= 1, nd = nd + 1; last = lev(t); end, end
    NJ = jL(max(Lv(:))) - jL(min(Lv(:)));
    dJles = mean(J(outL)) - mean(J(inL));
    fprintf(' %-4s %3d  %-6s %2d  %-5s %4.1f (%4.1f) |   %5.2f/%4.2f  %5.2f/%4.2f |   %5.3f/%5.3f  %5.3f/%5.3f  |  %4d    %4d   %5.0f        %5.1f\n', ...
            cmp, DR, mp, nb, mon, amb, ambc, cnr(P, inC, outC), cnr(P, inL, outL), cnr(J, inC, outC), cnr(J, inL, outL), ...
            gc(P, inC, outC), gc(P, inL, outL), gc(J, inC, outC), gc(J, inL, outL), numel(unique(P(:))), nd, NJ, dJles);
    if q == 1, Pref = P; end
    end

%% figure: reference case image, histograms per stage
figure('Position', [60 60 1200 420]);
subplot(1, 3, 1); imagesc(xc*1e3, zc*1e3, Pref, [0 255]); colormap(gray(256)); axis image;
xlabel('x (mm)'); ylabel('depth (mm)'); title('log 60 dB, linear map, 8 bit');
subplot(1, 3, 2);
Dd = 20*log10(A(tissue)/10^(D50/20));
hist(Dd(:), 80); xlabel('envelope (dB re tissue median)'); ylabel('pixels'); title('tissue envelope, dB');
subplot(1, 3, 3);
hist(Pref(tissue), 0:255); xlim([0 255]); xlabel('P-value'); title('tissue P-values (8 bit)');
