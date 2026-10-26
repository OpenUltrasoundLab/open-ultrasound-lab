%% m25_mv_esbmv.m  -  minimum variance (Capon) and eigenspace-based MV for ultrasound
% Purpose : (A) narrowband focal-plane analysis: MV vs DAS beam pattern with an
%           interferer, signal cancellation of a coherent echo versus subarray length L,
%           and the lone-point width law 0.914/sqrt(N*SNR); (B) broadband plane-wave
%           images of a point pair and a cyst in speckle with DAS, MV (L = 16, 32;
%           temporal averaging K = 0, 2) and ESBMV, with FWHM, pair dip, CR, gCNR,
%           speckle SNR and run time.
% Change  : dx_int, INR_dB (part A); L, Kt, Delta, delta, noise_dB (part B).
% Expect  : coherent echo cancelled (about -66 dB) at L = N, recovered (-0.3 dB) at
%           L = 48; MV point FWHM ~0.03 mm vs 0.48 mm DAS; K = 0 ruins speckle SNR.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8. Runtime ~10 s.
clear; close all; rng(2);
c = 1540; f0 = 5e6; lam = c/f0; N = 64; pitch = 0.3e-3; zF = 30e-3;
psi = @(dx) 2*pi*pitch*dx/(lam*zF);                  % phase step per element for a lateral offset dx at zF

%% A1. Beam pattern of MV weights (exact covariance, incoherent interferer)
dx_int = 2e-3; INR_dB = 10; SNR_dB = 20; L = 32; Delta = 1e-2; P = N - L + 1;
v = @(dx, n) exp(1j*psi(dx)*(0:n-1)');
R = ones(L) + 10^(INR_dB/10)*(v(dx_int,L)*v(dx_int,L)') + 10^(-SNR_dB/10)*eye(L);
R = R + Delta*real(trace(R))/L*eye(L); w = R\ones(L,1); w = w/sum(w);
xs = linspace(-8e-3, 8e-3, 1601); Bmv = zeros(size(xs)); Bdas = Bmv;
for k = 1:numel(xs)
  Bmv(k) = abs(w'*v(xs(k),L))*abs(mean(exp(1j*psi(xs(k))*(0:P-1))));   % subarray weights x Dirichlet of the smoothing
  Bdas(k) = abs(sum(v(xs(k),N)))/N;
end
[~, ki] = min(abs(xs - dx_int));
fprintf('A1 interferer at %.1f mm, %d dB: MV response %.1f dB, DAS %.1f dB (exact covariance)\n', dx_int*1e3, INR_dB, 20*log10(Bmv(ki)), 20*log10(Bdas(ki)));

%% A2. Signal cancellation: coherent interferer (same pulse), desired-signal gain vs L
dxc = 1.5e-3; beta = 1; u = 1 + beta*exp(1j*psi(dxc)*(0:N-1)'); Ls = 1:N; gain = zeros(size(Ls));
for iL = Ls
  P = N - iL + 1; R = zeros(iL);
  for l = 1:P, uu = u(l:l+iL-1); R = R + uu*uu'; end
  R = R/P + 10^(-SNR_dB/10)*eye(iL); R = R + Delta*real(trace(R))/iL*eye(iL);
  w = R\ones(iL,1); w = w/sum(w);
  gain(iL) = abs(1 + beta*(w'*exp(1j*psi(dxc)*(0:iL-1)'))*mean(exp(1j*psi(dxc)*(0:P-1))));
end
fprintf('A2 coherent echo at %.1f mm: desired-signal gain %.1f dB (L = 64), %.1f dB (L = 48), %.1f dB (L = 32), %.1f dB (L = 1)\n', dxc*1e3, 20*log10(gain([64 48 32 1])));

%% A3. Lone point: MV -6 dB width vs SNR (closed form vs direct solve)
xx = linspace(-0.3e-3, 0.3e-3, 6001); das = abs(sum(exp(1j*psi(linspace(-1e-3,1e-3,20001))'*(0:N-1)), 2))/N;
xd = linspace(-1e-3, 1e-3, 20001); wd = max(xd(das >= 0.5)) - min(xd(das >= 0.5));
for snr = [10 20 30]
  s2 = 10^(-snr/10); va = abs(sum(exp(1j*psi(xx)'*(0:N-1)), 2)); y = s2*va./(N*(s2 + N) - va.^2);
  Ri = inv(ones(N) + s2*eye(N)); vt = exp(1j*psi(0.05e-3)*(0:N-1)'); wt = Ri*vt/(vt'*Ri*vt); yd = abs(wt'*ones(N,1));
  [~, k5] = min(abs(xx - 0.05e-3)); ii = find(y >= 0.5); i1 = ii(1); i2 = ii(end);  % interpolate the -6 dB crossings (grid step 0.1 um)
  w6 = interp1(y([i2 i2+1]), xx([i2 i2+1]), 0.5) - interp1(y([i1-1 i1]), xx([i1-1 i1]), 0.5);
  fprintf('A3 SNR %2d dB: MV -6 dB width %.4f mm = %.4f x DAS (law 0.914/sqrt(N SNR) = %.4f); closed form vs solve at 0.05 mm: %.2e\n', snr, w6*1e3, w6/wd, 0.914/sqrt(N*10^(snr/10)), abs(y(k5) - yd));
end

%% B. Broadband plane-wave data: pair at z = 24 mm, cyst in speckle, channel noise
fs = 20e6; bw = 0.6; sigf = bw*f0/(2*sqrt(2*log(2))); sigt = 1/(2*pi*sigf); ang = (-4:4)*pi/180; Na = numel(ang);
xe = ((0:N-1) - (N-1)/2)*pitch; t0 = 25e-6; Nt = 560; tt = t0 + (0:Nt-1)'/fs; noise_dB = -20;
area = [-5.5 5.5 26 34]*1e-3; nsc = round(11*8*80); xsc = area(1) + 11e-3*rand(nsc,1); zsc = area(3) + 8e-3*rand(nsc,1);
keep = (xsc.^2 + (zsc-30e-3).^2) > (2.5e-3)^2; xsc = [xsc(keep); -0.25e-3; 0.25e-3; -2.5e-3]; zsc = [zsc(keep); 24e-3; 24e-3; 24e-3];
asc = [randn(sum(keep),1); 15; 15; 15];
IQ = zeros(Nt, N, Na); K4 = ceil(4*sigt*fs); E = repmat(1:N, numel(xsc), 1); tic;
for a = 1:Na
  tau = (xsc*sin(ang(a)) + zsc*cos(ang(a)))/c + sqrt((xsc - xe).^2 + zsc.^2)/c; ph = asc.*exp(-1j*2*pi*f0*tau); n0 = (tau - t0)*fs;
  accr = zeros(Nt*N,1); acci = accr;
  for k = -K4:K4
    n = floor(n0) + k; g = exp(-((n - n0)/fs).^2/(2*sigt^2)); ok = n >= 0 & n < Nt; idx = n(ok) + 1 + (E(ok)-1)*Nt; vv = g(ok).*ph(ok);
    accr = accr + accumarray(idx, real(vv), [Nt*N 1]); acci = acci + accumarray(idx, imag(vv), [Nt*N 1]);
  end
  IQ(:,:,a) = reshape(accr + 1j*acci, Nt, N);
end
Pref = mean(mean(abs(IQ(tt > 2*27.5e-3/c & tt < 2*33e-3/c, :, 5)).^2)); taps = exp(-((-K4:K4)'/fs).^2/(2*sigt^2));
for a = 1:Na, IQ(:,:,a) = IQ(:,:,a) + sqrt(Pref*10^(noise_dB/10)/sum(taps.^2))*conv2((randn(Nt+2*K4,N) + 1j*randn(Nt+2*K4,N))/sqrt(2), taps, 'valid'); end
fprintf('B  simulated %d scatterers in %.1f s\n', numel(xsc), toc);
dz = 0.075e-3; Kt = 2; delta = 0.5; Delta = 1e-2;
grids = {(-4:0.1:4)*1e-3, 23e-3 + (0:133)*dz; [(-0.8:0.005:0.8) (-2.8:0.005:-2.2)]*1e-3, 24e-3 + (-Kt:Kt)*dz};
for gi = 1:2
  xg = grids{gi,1}; zg = grids{gi,2}; [Zp, Xp] = ndgrid(zg, xg); Xp = Xp(:); Zp = Zp(:); S = zeros(numel(Xp), N);
  Ptx = exp(1j*2*pi*f0*(Xp*sin(ang) + Zp*cos(ang))/c); Ttx = (Xp*sin(ang) + Zp*cos(ang))/c;
  for e = 1:N
    trx = sqrt((Xp - xe(e)).^2 + Zp.^2)/c; acc = 0;
    for a = 1:Na
      fn = (Ttx(:,a) + trx - t0)*fs; n = floor(fn); fr = fn - n; ok = n >= 0 & n < Nt-1; col = IQ(:,e,a); vv = zeros(size(fn));
      vv(ok) = col(n(ok)+1).*(1-fr(ok)) + col(n(ok)+2).*fr(ok); acc = acc + vv.*Ptx(:,a);
    end
    S(:,e) = acc.*exp(1j*2*pi*f0*trx);
  end
  S = reshape(S, numel(zg), numel(xg), N); nz = numel(zg); nx = numel(xg);
  if gi == 1, rows = 1:nz; else, rows = Kt+1; end
  out = zeros(nz, nx, 5); tm = zeros(1,5);                % DAS, MV L32 K2, MV L32 K0, MV L16 K2, ESBMV L32 K2
  out(:,:,1) = abs(sum(S,3))/N;
  cfg = [32 Kt; 32 0; 16 Kt];
  for ix = 1:nx
    for r = rows
      for ci = 1:3
        t1 = tic; L = cfg(ci,1); K = cfg(ci,2); P = N - L + 1; id = (0:L-1)' + (1:P); X = [];
        for k = -K:K, if r+k >= 1 && r+k <= nz, x = squeeze(S(r+k, ix, :)); X = [X, x(id)]; end, end %#ok<AGROW>
        Rm = (X*X')/size(X,2); w = (Rm + Delta*real(trace(Rm))/L*eye(L))\ones(L,1); w = w/sum(w);
        x = squeeze(S(r, ix, :)); b = mean(x(id), 2); out(r, ix, 1+ci) = abs(w'*b); tm(1+ci) = tm(1+ci) + toc(t1);
        if ci == 1
          t2 = tic; [V, D] = eig((Rm + Rm')/2); [lv, o] = sort(real(diag(D)), 'descend'); Es = V(:, o(lv >= delta*lv(1)));
          out(r, ix, 5) = abs((Es*(Es'*w))'*b); tm(5) = tm(5) + toc(t2);
        end
      end
    end
  end
  if gi == 1, img = out; xi = xg; zi = zg; T = tm; T(5) = T(5) + T(2); else, prof = squeeze(out(Kt+1, :, :)); xp = xg; end
end

%% Metrics
names = {'DAS', 'MV L32 K2', 'MV L32 K0', 'MV L16 K2', 'ESBMV'};
[ZI, XI] = ndgrid(zi, xi); cyst = (XI.^2 + (ZI-30e-3).^2) < (2e-3)^2; bgr = abs(XI) >= 3e-3 & abs(XI) <= 3.9e-3 & ZI >= 28e-3 & ZI <= 32e-3;
fprintf('\n%-10s %7s %7s %6s %6s %6s %7s\n', 'method', 'FWHM', 'dip', 'CR', 'gCNR', 'sSNR', 'time');
for m = 1:5
  A = img(:,:,m); ci = A(cyst); bi = A(bgr); v = [ci; bi]; lab = [zeros(numel(ci),1); ones(numel(bi),1)]; [~, o] = sort(v); lab = lab(o);
  bins = floor((0:numel(v)-1)'*100/numel(v)); g = 1 - sum(min(accumarray(bins(lab==0)+1,1,[100 1])/numel(ci), accumarray(bins(lab==1)+1,1,[100 1])/numel(bi)));
  pr = prof(:,m); sel = xp < -2e-3; q = pr(sel); xq = xp(sel); [pk, ip] = max(q); above = xq(q >= pk/2);
  pp = pr(~sel); xx = xp(~sel); p1 = max(pp(abs(xx + 0.25e-3) < 0.05e-3)); p2 = max(pp(abs(xx - 0.25e-3) < 0.05e-3));
  dip = 20*log10(min(pp(abs(xx) < 0.15e-3))/min(p1, p2));
  fprintf('%-10s %6.3fmm %6.1fdB %6.1f %6.3f %6.2f %6.1fs\n', names{m}, (max(above) - min(above))*1e3, dip, 10*log10(mean(bi.^2)/mean(ci.^2)), g, mean(bi)/std(bi,1), T(m));
end

%% Figures
figure('Position', [50 50 1200 420]);
subplot(1,2,1); plot(xs*1e3, 20*log10(Bdas), 'k', xs*1e3, 20*log10(Bmv + 1e-12), 'b'); hold on; plot([1 1]*dx_int*1e3, [-80 0], 'r--');
ylim([-80 2]); xlabel('lateral position at z = 30 mm (mm)'); ylabel('beam pattern (dB)'); legend('DAS', 'MV (L = 32)', 'interferer'); title('A1: MV places a null on the interferer');
subplot(1,2,2); plot(Ls, 20*log10(gain), 'b', 'LineWidth', 1.5); xlabel('subarray length L'); ylabel('desired-signal gain (dB)'); title('A2: coherent echo cancels without smoothing'); grid on;
figure('Position', [50 50 1300 330]);
for m = 1:5
  subplot(1,5,m); A = img(:,:,m)/max(max(img(:,:,m))); imagesc(xi*1e3, zi*1e3, 20*log10(A + 1e-12), [-60 0]); axis image; colormap(gray);
  title(names{m}); xlabel('x (mm)'); ylabel('z (mm)');
end
