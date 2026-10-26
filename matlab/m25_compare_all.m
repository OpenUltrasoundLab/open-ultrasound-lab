%% m25_compare_all.m  -  ten beamformers on one simulated plane-wave dataset
% Purpose : simulate 64-channel IQ data for 9 compounded plane waves (point targets,
%           a 0.5 mm pair, an anechoic cyst in speckle, a bright scatterer, channel
%           noise, optional receive phase screen), beamform with DAS, MV, ESBMV, CF,
%           GCF, PCF, Wiener, F-DMAS, SLSC and FXPF, and print a metric table:
%           FWHM, PSL, pair dip, CR, CNR and gCNR at two display ranges, speckle SNR, time.
% Change  : noise_dB (channel noise re speckle), screen_ns (rms phase screen),
%           L, Delta, K (MV), delta (ESBMV), M0 (GCF), M (SLSC), p, W (FXPF).
% Expect  : MV/ESBMV FWHM ~0.03 mm and a deep pair dip; CF/PCF/ESBMV CR far above
%           DAS but gCNR close to DAS; SLSC lowest CR, highest gCNR; runtime ~20 s.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

%% Acquisition and phantom
c = 1540; f0 = 5e6; fs = 20e6; lam = c/f0;            % m/s, Hz, Hz, m
N = 64; pitch = 0.3e-3; xe = ((0:N-1) - (N-1)/2)*pitch; % elements (m)
bw = 0.6; sigf = bw*f0/(2*sqrt(2*log(2))); sigt = 1/(2*pi*sigf);   % Gaussian pulse
ang = (-4:4)*pi/180; Na = numel(ang);                 % plane waves, 1 deg steps
t0 = 25e-6; Nt = 580; tt = t0 + (0:Nt-1)'/fs;         % IQ time axis (s)
noise_dB = -20; screen_ns = 0;                        % channel noise, phase screen
% speckle (80 /mm^2) with an anechoic cyst, three points, one bright scatterer
area = [-6.8 6.8 26 34.6]*1e-3; dens = 80e6;          % m, scatterers per m^2
nsc = round((area(2)-area(1))*(area(4)-area(3))*dens);
xs = area(1) + (area(2)-area(1))*rand(nsc,1); zs = area(3) + (area(4)-area(3))*rand(nsc,1);
keep = (xs.^2 + (zs-30e-3).^2) > (2.5e-3)^2; xs = xs(keep); zs = zs(keep); as = randn(numel(xs),1);
xs = [xs; -3e-3; 2.25e-3; 2.75e-3; 4.2e-3]; zs = [zs; 24e-3; 24e-3; 24e-3; 30.5e-3]; as = [as; 15; 15; 15; 60];

%% Channel data: iq(t) = a g(t - tau) exp(-j 2 pi f0 tau), accumulated sample by sample
tic; IQ = zeros(Nt, N, Na); K4 = ceil(4*sigt*fs); E = repmat(1:N, numel(xs), 1);
for a = 1:Na
  tau = (xs*sin(ang(a)) + zs*cos(ang(a)))/c + sqrt((xs - xe).^2 + zs.^2)/c;
  ph = as.*exp(-1j*2*pi*f0*tau); n0 = (tau - t0)*fs; accr = zeros(Nt*N,1); acci = accr;
  for k = -K4:K4
    n = floor(n0) + k; g = exp(-((n - n0)/fs).^2/(2*sigt^2)); ok = n >= 0 & n < Nt;
    idx = n(ok) + 1 + (E(ok)-1)*Nt; v = g(ok).*ph(ok);
    accr = accr + accumarray(idx, real(v), [Nt*N 1]); acci = acci + accumarray(idx, imag(v), [Nt*N 1]);
  end
  IQ(:,:,a) = reshape(accr + 1j*acci, Nt, N);
end
band = tt > 2*27.5e-3/c & tt < 2*33e-3/c; Pref = mean(mean(abs(IQ(band,:,5)).^2));
taps = exp(-((-K4:K4)'/fs).^2/(2*sigt^2)); sn = sqrt(Pref*10^(noise_dB/10)/sum(taps.^2));
for a = 1:Na
  w = (randn(Nt+2*K4, N) + 1j*randn(Nt+2*K4, N))/sqrt(2);
  IQ(:,:,a) = IQ(:,:,a) + sn*conv2(w, taps, 'valid');
end
scr = zeros(1,N);                                      % receive phase screen (s)
if screen_ns > 0
  wv = randn(1, N+40); ker = exp(-(-20:20).^2/(2*(1.5e-3/pitch)^2)); scr = conv(wv, ker, 'valid');
  scr = scr - mean(scr); scr = scr/sqrt(mean(scr.^2))*screen_ns*1e-9;
end
fprintf('simulated %d scatterers, 9 angles x 64 channels in %.1f s\n', numel(xs), toc);

%% Parameters
L = 32; Delta = 1e-2; Kt = 2; delta = 0.5; M0 = 1; gam = 1; M = 8; Ks = 8; p = 4; W = 16; mu = 1e-3;
sub = 4; dz = 0.075e-3; dzf = dz/sub; fsz = c/(2*dzf);  % output and RF axial steps, RF rate (Hz)
names = {'DAS','MV','ESBMV','CF','GCF','PCF','Wiener','F-DMAS','SLSC','FXPF'}; nm = numel(names);
T = zeros(1, nm);                                     % accumulated runtime per method (s)

%% Two passes: pass 1 = image grid, pass 2 = lateral profile through the point at (-3, 24) mm
for pass = 1:2
  if pass == 1
    xg = (-6:0.15:6)*1e-3; z0 = 22e-3; nzo = 161; rows = 1:nzo;
  else
    xg = [(-5.48:0.04:-3.64) (-3.6:0.02:-3.16) (-3.15:0.004:-2.85) (-2.84:0.02:-2.4) (-2.36:0.04:-0.52) (1.5:0.02:3.5)]*1e-3;
    z0 = 22.8e-3; nzo = 33; rows = 17;
  end
  zf = z0 + (0:(nzo-1)*sub)'*dzf; nzf = numel(zf); nx = numel(xg);
  % aperture data: S(row, col, element), delayed, phase-corrected, compounded (linear interpolation of IQ)
  [Zp, Xp] = ndgrid(zf, xg); Xp = Xp(:); Zp = Zp(:); np_ = numel(Xp); S = zeros(np_, N);
  Ttx = (Xp*sin(ang) + Zp*cos(ang))/c; Ptx = exp(1j*2*pi*f0*Ttx);        % np x Na
  for e = 1:N
    trx = sqrt((Xp - xe(e)).^2 + Zp.^2)/c; acc = zeros(np_, 1);
    for a = 1:Na
      fn = (Ttx(:,a) + trx - scr(e) - t0)*fs; n = floor(fn); fr = fn - n; ok = n >= 0 & n < Nt-1;
      col = IQ(:,e,a); v = zeros(np_, 1); v(ok) = col(n(ok)+1).*(1 - fr(ok)) + col(n(ok)+2).*fr(ok);
      acc = acc + v.*Ptx(:,a);
    end
    S(:,e) = acc.*exp(1j*2*pi*f0*(trx - scr(e)));
  end
  S = reshape(S, nzf, nx, N); So = S(1:sub:end,:,:);  % output-row samples
  out = zeros(nzo, nx, nm);
  % --- DAS and coherence weights (per sample)
  tic; ys = sum(So,3); out(:,:,1) = abs(ys)/N; T(1) = T(1) + toc;
  tic; P = sum(abs(So).^2,3); CF = abs(ys).^2./(N*P + eps); out(:,:,4) = CF.*abs(ys)/N; T(4) = T(4) + toc;
  tic; Sk = fft(So, [], 3); low = abs(Sk(:,:,1)).^2;
  for k = 1:M0, low = low + abs(Sk(:,:,1+k)).^2 + abs(Sk(:,:,N+1-k)).^2; end
  out(:,:,5) = low./(N*P + eps).*abs(ys)/N; T(5) = T(5) + toc;
  tic; phs = angle(So); pha = phs - pi*sign(phs); sg = min(std(phs,1,3), std(pha,1,3));
  out(:,:,6) = max(0, 1 - gam*sg/(pi/sqrt(3))).*abs(ys)/N; T(6) = T(6) + toc;
  tic; out(:,:,7) = (N*CF./((N-1)*CF + 1)).*abs(ys)/N; T(7) = T(7) + toc;
  % --- F-DMAS on RF = real part, fine axial grid, Gaussian band-pass at 2 f0
  tic; R = real(S); sh = sign(R).*sqrt(abs(R)); yd = (sum(sh,3).^2 - sum(abs(R),3))/2;
  nf = 2^nextpow2(nzf + 64); f = [0:nf/2, -nf/2+1:-1]'*fsz/nf; H = 2*exp(-(f - 2*f0).^2/(2*(sqrt(2)*sigf)^2)).*(f > 0);
  env = abs(ifft(fft(yd, nf, 1).*H, [], 1)); out(:,:,8) = env(1:sub:nzf, :)/(N*(N-1)/2); T(8) = T(8) + toc;
  % --- SLSC on RF, axial kernel +-Ks fine samples (kernel sums by cumulative sums along depth)
  tic; ro = (rows-1)*sub + 1; hi = min(nzf, ro + Ks); lo = max(1, ro - Ks);
  Cs = cumsum(cat(1, zeros(1, nx, N), R.^2), 1); En = Cs(hi+1,:,:) - Cs(lo,:,:); sl = zeros(numel(rows), nx);
  for m = 1:M
    Cs = cumsum(cat(1, zeros(1, nx, N-m), R(:,:,1:N-m).*R(:,:,1+m:N)), 1); num = Cs(hi+1,:,:) - Cs(lo,:,:);
    sl = sl + mean(num./sqrt(En(:,:,1:N-m).*En(:,:,1+m:N) + eps), 3);
  end
  out(rows, :, 9) = max(sl, 0); T(9) = T(9) + toc;
  % --- MV and ESBMV: subarrays L, temporal averaging +-Kt output samples, loading Delta
  Pn = N - L + 1; idx = (0:L-1)' + (1:Pn); tmv = 0; tes = 0;
  for ix = 1:nx
    for r = rows
      t1 = tic; X = [];
      for k = -Kt:Kt
        if r+k >= 1 && r+k <= nzo, x = squeeze(So(r+k, ix, :)); X = [X, x(idx)]; end %#ok<AGROW>
      end
      Rm = (X*X')/size(X,2); Rl = Rm + Delta*real(trace(Rm))/L*eye(L);
      w = Rl\ones(L,1); w = w/(ones(1,L)*w); x = squeeze(So(r, ix, :)); b = mean(x(idx), 2);
      out(r, ix, 2) = abs(w'*b); tmv = tmv + toc(t1);
      t2 = tic; [V, D] = eig((Rm + Rm')/2); lam_ = real(diag(D)); [lam_, o] = sort(lam_, 'descend'); V = V(:,o);
      Es = V(:, lam_ >= delta*lam_(1)); we = Es*(Es'*w); out(r, ix, 3) = abs(we'*b); tes = tes + toc(t2);
    end
  end
  T(2) = T(2) + tmv; T(3) = T(3) + tmv + tes;
  % --- FXPF: Hann windows (W samples, hop W/2); per bin, forward and backward AR(p) least squares across the
  %     aperture, solved for all columns at once (vectorized elimination); predictions averaged; then DAS
  tic; win = 0.5 - 0.5*cos(2*pi*(0:W-1)'/W); hop = W/2; Y = zeros(nzo, nx, N);
  for s0 = -hop:hop:nzo-1
    rr = s0 + (1:W); ok = rr >= 1 & rr <= nzo; seg = zeros(W, nx, N); seg(ok,:,:) = So(rr(ok),:,:);
    F = fft(seg.*win, [], 1); G = zeros(W, nx, N);
    for kb = 1:W
      Xb = reshape(F(kb,:,:), nx, N); pr = zeros(nx, N); cnt = zeros(1, N);
      for dirn = [1 -1]
        if dirn > 0, tgt = p+1:N; else, tgt = 1:N-p; end
        A = zeros(nx, numel(tgt), p); for q = 1:p, A(:,:,q) = Xb(:, tgt - dirn*q); end
        Gm = zeros(nx, p, p); h = zeros(nx, p);
        for q1 = 1:p, h(:,q1) = sum(conj(A(:,:,q1)).*Xb(:,tgt), 2); for q2 = 1:p, Gm(:,q1,q2) = sum(conj(A(:,:,q1)).*A(:,:,q2), 2); end, end
        tr = zeros(nx,1); for q = 1:p, tr = tr + real(Gm(:,q,q)); end
        for q = 1:p, Gm(:,q,q) = Gm(:,q,q) + mu*tr/p + 1e-30; end
        for k = 1:p                                   % Gaussian elimination, vectorized over columns
          for i = k+1:p
            fct = Gm(:,i,k)./Gm(:,k,k); for j = k:p, Gm(:,i,j) = Gm(:,i,j) - fct.*Gm(:,k,j); end, h(:,i) = h(:,i) - fct.*h(:,k);
          end
        end
        cf_ = zeros(nx, p);
        for k = p:-1:1, sres = h(:,k); for j = k+1:p, sres = sres - Gm(:,k,j).*cf_(:,j); end, cf_(:,k) = sres./Gm(:,k,k); end
        P_ = zeros(nx, numel(tgt)); for q = 1:p, P_ = P_ + cf_(:,q).*A(:,:,q); end
        pr(:, tgt) = pr(:, tgt) + P_; cnt(tgt) = cnt(tgt) + 1;
      end
      G(kb,:,:) = reshape(pr./cnt, [1 nx N]);
    end
    g = ifft(G, [], 1); Y(rr(ok),:,:) = Y(rr(ok),:,:) + g(ok,:,:);
  end
  out(:,:,10) = abs(sum(Y,3))/N; T(10) = T(10) + toc;
  if pass == 1, img = out; xi = xg; zi = z0 + (0:nzo-1)*dz; else, prof = squeeze(out(17,:,:)); xp = xg; end
  fprintf('pass %d: %d x %d pixels, cumulative times (s): %s\n', pass, nzo, nx, sprintf('%.1f ', T));
end

%% Metrics
[ZI, XI] = ndgrid(zi, xi); cyst = (XI.^2 + (ZI-30e-3).^2) < (2e-3)^2;
bgr = XI >= -5.8e-3 & XI <= -3.4e-3 & ZI >= 28e-3 & ZI <= 32e-3;
fprintf('\n%-7s %6s %6s %6s %6s %7s %7s %6s %7s %7s %6s %7s\n', 'method', 'FWHM', 'PSL', 'dip', 'CR', 'CNR60', 'CNR40', 'gCNR', 'gCNR60', 'gCNR40', 'sSNR', 'time');
for m = 1:nm
  A = max(img(:,:,m), 0); mx = max(A(:)); ci = max(A(cyst), 1e-12*mx); bi = max(A(bgr), 1e-12*mx);
  res = zeros(1, 5); DRs = [Inf 60 40];
  for d = 1:3
    vc = max(20*log10(ci/mx), -DRs(d)); vb = max(20*log10(bi/mx), -DRs(d));
    v = [vc; vb]; lab = [zeros(numel(vc),1); ones(numel(vb),1)]; [vs, o] = sort(v); lab = lab(o); n = numel(v);
    bins = floor((0:n-1)'*100/n); [~, first, ic] = unique(vs, 'first'); bins = bins(first(ic));
    hc = accumarray(bins(lab == 0)+1, 1, [100 1])/numel(vc); hb = accumarray(bins(lab == 1)+1, 1, [100 1])/numel(vb);
    res(d) = 1 - sum(min(hc, hb));
    if d > 1, res(3+d-1) = abs(mean(vb) - mean(vc))/sqrt(var(vb,1) + var(vc,1)); end
  end
  CR = 10*log10(mean(bi.^2)/mean(ci.^2)); sSNR = mean(bi)/std(bi,1);
  sel = xp < 0; pr = prof(sel, m); xq = xp(sel); [pk, ip] = max(pr.*(abs(xq(:) + 3e-3) < 0.6e-3)); half = pk/2;
  il = ip; while il > 1 && pr(il) > half, il = il - 1; end
  ir = ip; while ir < numel(pr) && pr(ir) > half, ir = ir + 1; end
  xl = xq(il) + (half - pr(il))/(pr(il+1) - pr(il))*(xq(il+1) - xq(il));
  xr = xq(ir-1) + (half - pr(ir-1))/(pr(ir) - pr(ir-1))*(xq(ir) - xq(ir-1));
  ml = il; while ml > 1 && pr(ml-1) <= pr(ml), ml = ml - 1; end
  mr = ir; while mr < numel(pr) && pr(mr+1) <= pr(mr), mr = mr + 1; end
  side = pr; side(ml:mr) = 0; side(abs(xq(:) - xq(ip)) >= 2.5e-3) = 0; PSL = 20*log10(max(max(side), 1e-9*pk)/pk);
  q = prof(~sel, m); xx = xp(~sel); p1 = max(q(abs(xx - 2.25e-3) < 0.1e-3)); p2 = max(q(abs(xx - 2.75e-3) < 0.1e-3));
  dip = 20*log10(min(q(xx > 2.3e-3 & xx < 2.7e-3))/min(p1, p2));
  fprintf('%-7s %6.3f %6.1f %6.1f %6.1f %7.2f %7.2f %6.3f %7.3f %7.3f %6.2f %6.1fs\n', names{m}, (xr-xl)*1e3, PSL, dip, CR, res(4), res(5), res(1), res(2), res(3), sSNR, T(m));
end

%% Figure: 2 x 5 image grid, 60 dB (SLSC linear)
figure('Position', [50 50 1300 560]);
for m = 1:nm
  subplot(2, 5, m); A = max(img(:,:,m), 0); A = A/max(A(:));
  if m == 9, imagesc(xi*1e3, zi*1e3, A, [0 1]); else, imagesc(xi*1e3, zi*1e3, 20*log10(A + 1e-12), [-60 0]); end
  axis image; colormap(gray); title(names{m}); xlabel('x (mm)'); ylabel('z (mm)');
end
