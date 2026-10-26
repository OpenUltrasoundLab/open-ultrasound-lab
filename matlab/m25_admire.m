%% m25_admire.m  --  ADMIRE (aperture domain model image reconstruction) on one image line
% Purpose : a compact version of the ADMIRE pipeline of Module 25 (Byram et al., IEEE TUFFC 2015;
%           model space and regularization from the authors' released code, VU-BEAM-Lab/ADMIRE):
%           delayed aperture data -> short-time FFT along depth -> per window and frequency, fit the
%           aperture vector with modelled wavefronts of scatterers inside (ROI) and outside (off-axis)
%           the resolution cell (elastic net) -> keep the ROI part -> inverse FFT -> DAS.
% Scene   : one line at x = 0 through a point target at (0, 30) mm, 120 off-axis scatterers at
%           1.5 <= |x| <= 9 mm (z = 28-32 mm) and white channel noise; 64 elements, 0.3 mm pitch,
%           5 MHz, 60 % bandwidth, one 0 deg plane wave, dynamic receive focusing.
% Differences from the page: one plane wave instead of nine, the reduced model space of the page's
%           live button (about 900 predictors), and the elastic net solved by FISTA (vectorized)
%           instead of coordinate descent; the objective, alpha and lambda rule are the same.
% Change  : clutAmp (clutter strength), snrDB, cLam (lambda scale), useOuter (0 = no off-axis model).
% Expected: printed target loss (ADMIRE vs DAS) and clutter level relative to the target, both
%           methods; figure with the depth profiles and the aperture fit at the target window.
% Runtime : about 10-20 s in Octave 8.
rng(1);
c = 1.54; f0 = 5; BW = 0.6; N = 64; pitch = 0.3;            % mm/us, MHz, -, -, mm
xe = ((0:N-1) - (N-1)/2) * pitch;                           % element x (mm)
sigf = BW*f0/(2*sqrt(2*log(2))); sigt = 1/(2*pi*sigf);      % pulse: Gaussian envelope (us)
dz = 0.075; z = (27.5:dz:32.5)'; nz = numel(z);             % image rows on the line x = 0 (mm)
clutAmp = 0.25; nClut = 120; snrDB = 20; cLam = 0.0189; alpha = 0.9; useOuter = 1;

%% Delayed aperture data on the line x = 0: s(n,e) = sum_s a_s g(dt) exp(j 2 pi f0 dt), dt = tau_pix - tau_s
xc = (1.5 + 7.5*rand(nClut,1)) .* sign(rand(nClut,1) - 0.5); zcl = 28 + 4*rand(nClut,1);
ac = clutAmp * (randn(nClut,1) + 1i*randn(nClut,1)) / sqrt(2);
tauPix = (z + sqrt(xe.^2 + z.^2)) / c;                      % nz x N, 0 deg transmit + receive
mk = @(xs, zs, as) reshape(sum(reshape(as,1,1,[]) .* exp(-(tauPix - reshape((zs + sqrt((xs - xe).^2 + zs.^2))', 1, N, []) / c).^2 ...
     / (2*sigt^2)) .* exp(1i*2*pi*f0*(tauPix - reshape((zs + sqrt((xs - xe).^2 + zs.^2))', 1, N, []) / c)), 3), nz, N);
sT = mk(0, 30, 1);                                          % target only
sC = zeros(nz, N); for k = 1:30:nClut, idx = k:min(k+29, nClut); sC = sC + mk(xc(idx), zcl(idx), ac(idx)); end
sigP = mean(abs(sT(abs(z-30) < 0.15, :)).^2, 1); nsd = sqrt(mean(sigP) * 10^(-snrDB/10) / 2);
sN = nsd * (randn(nz, N) + 1i*randn(nz, N));               % white noise, SNR re target peak per channel
scenes = {sT + sC + sN, sC + sN, sT};                       % A: target + clutter + noise, B: clutter + noise, C: target

%% STFT and model parameters (released rule: window = ceil(FWHM * fs), zero-pad x2, bins in f0 +- 0.6 BW f0)
fsa = c/(2*dz); Lw = ceil(sqrt(8*log(2))*sigt*fsa); Lp = 2*Lw;
fk = (0:Lp-1)*fsa/Lp; bins = find(fk >= f0*(1-0.6*BW) & fk <= f0*(1+0.6*BW));
kcal = [1.10 0.95 0.88];                                    % per-bin wavenumber calibration (page, fitted to this pulse)
h = Lw/fsa*c/2; st = c*sigt; nw = floor(nz/Lw);
gridv = @(a, s, b) a + (0:floor((b - a)/s + 1e-9))*s;      % MATLAB a:s:b without round-off surprises
outD = cell(3,1); fitShow = [];
tic;
for w = 1:nw
  rows = (w-1)*Lw + (1:Lw); zf = z(rows(1)); zc = mean(z(rows)); z2a = 2*zf; z2b = 2*z(rows(end));
  S = cell(3,1); for q = 1:3, X = fft(scenes{q}(rows, :), Lp, 1); S{q} = X; end
  R = cell(3,1); for q = 1:3, R{q} = zeros(Lp, N); end
  for b = 1:numel(bins)
    f = fk(bins(b)); kw = 2*pi*f*kcal(b)/c; lam = 2*pi/kw;
    resLat = 1.2067/(N*pitch)*lam*zc; resAx = 2*resLat;    % FWHM of |sinc| aperture transform (release computes it by FFT)
    % predictor grid: ROI and outer (off-axis) scatterers, offsets d (reduced sampling of the page's live model)
    [xr, zr, dr] = ndgrid(gridv(-0.5, 0.1432, 0.5)*resLat, zc + gridv(-0.5, 0.572, 0.5)*resAx, gridv(-0.8, 0.388*lam, 0.4));
    L = pitch*N/2 + 1; [xq, zq, dq] = ndgrid(gridv(-L, 1.4228*resLat, L), zc + gridv(-0.5, 0.7114, 0.5)*resAx, gridv(-8, 0.9688*lam, 3.2));
    keepO = (xq(:)/(resLat + 0.5)).^2 + ((zq(:) - zc)/(resAx + 0.5)).^2 > 1;
    xp = [xr(:); xq(keepO)]; zp = [zr(:); zq(keepO)]; dp = [dr(:); dq(keepO)]; nRoi = numel(xr);
    if ~useOuter, xp = xp(1:nRoi); zp = zp(1:nRoi); dp = dp(1:nRoi); end
    d0 = sqrt((xe - xp).^2 + zp.^2) + dp + 2*zc - zp;     % P x N total path
    zd = (d0.^2 - xe.^2) ./ d0;                             % two-way position after dynamic receive focusing
    pmin = max(zd - h, z2a); pmax = min(zd + h, z2b);
    g = erf((zd - pmin)/st) - erf((zd - pmax)/st); g(zd - h > z2b | zd + h < z2a) = 0;
    ok = max(g, [], 2) > 0; g = g(ok, :); zd = zd(ok, :); isRoi = (1:numel(ok))' <= nRoi; isRoi = isRoi(ok);
    A = sqrt(g ./ max(g, [], 2)) .* exp(-1i*kw*zd);        % modelled aperture signals (rows)
    A = (A ./ sqrt(sum(abs(A).^2, 2))).';                   % N x P, unit-norm columns
    Xr = [real(A) -imag(A); imag(A) real(A)]; P = size(A, 2); No = 2*N;
    Lip = norm(Xr)^2/No;                                    % Lipschitz constant of the smooth part
    for q = 1:3
      y = S{q}(bins(b), :).'; yv = [real(y); imag(y)]; sd = std(yv, 1);
      if sd == 0, continue; end
      ys = yv/sd; lm = cLam*sqrt(mean(yv.^2))/sd; tstep = 1/(Lip + lm*(1-alpha));
      B = zeros(2*P, 1); V = B; tk = 1;                     % FISTA on (1/2No)||y - Xb||^2 + lm(alpha|b|_1 + (1-alpha)/2 |b|^2)
      for it = 1:400
        G = V - tstep*(Xr'*(Xr*V - ys)/No + lm*(1-alpha)*V);
        Bn = sign(G).*max(abs(G) - tstep*lm*alpha, 0);
        tn = (1 + sqrt(1 + 4*tk^2))/2; V = Bn + (tk - 1)/tn*(Bn - B);
        if norm(Bn - B) < 1e-7*max(1, norm(B)), B = Bn; break; end
        B = Bn; tk = tn;
      end
      B = B*sd; bc = B(1:P) + 1i*B(P+1:end);
      R{q}(bins(b), :) = (A(:, isRoi) * bc(isRoi)).';      % keep only the in-ROI part
      if q == 1 && abs(zc - 30) < Lw*dz/2 && b == 2, fitShow = struct('y', y, 'roi', A(:, isRoi)*bc(isRoi), 'out', A(:, ~isRoi)*bc(~isRoi), 'P', P, 'nRoi', sum(isRoi)); end
    end
  end
  for q = 1:3, r = ifft(R{q}, [], 1); outD{q}(rows, 1) = abs(sum(r(1:Lw, :), 2))/N; end
end
tRun = toc;
das = cellfun(@(s) abs(sum(s, 2))/N, scenes, 'UniformOutput', false);
nr = nw*Lw; zz = z(1:nr);

%% Numbers: target kept? clutter removed?
iT = abs(zz - 30) < 0.04;                                   % rows at the target depth
iC = abs(zz - 30) > 0.6;                                    % rows away from the target (clutter only)
pk = @(v) max(v(iT));
dasT = pk(das{3}(1:nr)); admT = pk(outD{3});               % target alone
fprintf('Window %d samples (%.2f mm), bins %s MHz, last fit %d predictors (%d ROI), %.1f s\n', Lw, Lw*dz, mat2str(round(fk(bins)*100)/100), fitShow.P, fitShow.nRoi, tRun);
fprintf('Target only: ADMIRE/DAS peak = %.2f dB (amplitude loss from the elastic-net shrinkage and the 3-bin band)\n', 20*log10(admT/dasT));
fprintf('Target + clutter + noise at the target: DAS %.2f dB, ADMIRE %.2f dB re its own target-only peak\n', 20*log10(pk(das{1}(1:nr))/dasT), 20*log10(pk(outD{1})/admT));
cD = 10*log10(mean(das{2}(iC).^2)/dasT^2); pA = mean(outD{2}(iC).^2);
if pA > 0, sA = sprintf('%.1f dB', 10*log10(pA/admT^2)); else, sA = 'exactly 0 (no ROI coefficient survives)'; end
fprintf('Clutter + noise only, rows > 0.6 mm from the target: DAS %.1f dB re target peak, ADMIRE %s\n', cD, sA);
fprintf('Clutter-only rows set exactly to zero by ADMIRE: %.0f %%; in the target window it keeps %.0f %% of the rows nonzero\n', 100*mean(outD{2}(iC) == 0), 100*mean(outD{1}(iT) > 0));

%% Figure
figure(1); clf;
subplot(1,2,1);
plot(zz, 20*log10(das{1}(1:nr)/dasT + 1e-6), 'k', zz, 20*log10(outD{1}/admT + 1e-6), 'r', ...
     zz, 20*log10(das{2}(1:nr)/dasT + 1e-6), 'k:', zz, 20*log10(outD{2}/admT + 1e-6), 'r:', 'LineWidth', 1.2);
ylim([-70 5]); xlabel('depth z (mm)'); ylabel('amplitude (dB re target-only peak)');
legend('DAS, target+clutter', 'ADMIRE, target+clutter', 'DAS, clutter only', 'ADMIRE, clutter only', 'Location', 'southwest');
title('Line x = 0, depth profiles');
subplot(1,2,2);
plot(xe, real(fitShow.y), 'k.', xe, real(fitShow.roi), 'r', xe, real(fitShow.out), 'b', 'LineWidth', 1.2);
xlabel('element x (mm)'); ylabel('Re Y(f) at 5.13 MHz (a.u.)');
legend('data', 'ROI part (kept)', 'outer part (removed)', 'Location', 'best');
title('Aperture fit, target window');
