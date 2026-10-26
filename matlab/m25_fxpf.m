%% m25_fxpf.m  -  f-x prediction filtering (FXPF) of delayed channel data
% Purpose : (1) show that q plane waves across the aperture obey an order-q linear
%           recursion (exact prediction for p >= q); (2) Monte Carlo of the FXPF gain on
%           white noise and on a focused echo versus channel SNR (regression dilution);
%           (3) image a cyst/point phantom with DAS and with FXPF (Hann windows of W
%           samples, 50 % overlap, per-bin forward + backward AR(p) least squares, DAS).
% Change  : p, W, mu, noise_dB; add more plane waves in section 1.
% Expect  : section 1 residual ~1e-15 once p >= q; noise-only gain about -14 dB for p = 4
%           (forward-only projection estimate 10 log10(p/(N-p)) = -11.8 dB); a 0 dB echo
%           loses ~2 dB (p = 4) and ~6 dB (p = 1); at -5 dB channel noise FXPF raises CR
%           and gCNR while the point width and speckle SNR stay. Runtime ~15 s in Octave.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(3);
N = 64; c = 1540; f0 = 5e6; lam = c/f0; pitch = 0.3e-3; z = 30e-3; mu = 1e-3;
psi = @(dx) 2*pi*pitch*dx/(lam*z);                      % phase step per element for a wave from lateral offset dx

%% 1. Linear prediction across the aperture: residual of a forward AR(p) fit to q plane waves
fprintf('relative LS residual of a forward AR(p) fit, noise-free, q plane waves (offsets 0, 1.3, -2.1 mm)\n');
offs = [0 1.3e-3 -2.1e-3]; amps = [1 0.7 0.5];
for q = 1:3
  x = zeros(N,1); for k = 1:q, x = x + amps(k)*exp(1j*(psi(offs(k))*(0:N-1)' + k)); end
  r = zeros(1,4);
  for p = 1:4
    A = zeros(N-p, p); for k = 1:p, A(:,k) = x(p+1-k:N-k); end; b = x(p+1:N); a = A\b; r(p) = norm(b - A*a)/norm(b);
  end
  fprintf('  q = %d: p = 1..4 -> %s\n', q, sprintf('%9.1e', r));
end

%% 2. Monte Carlo: FXPF prediction (forward + backward, loading mu) on a constant echo plus white noise
p_list = [1 2 4 8]; snr_list = -20:5:30; T = 200; Tn = 2000;   % trials per echo SNR, trials for noise only (heavy-tailed: needs more)
gs = zeros(numel(p_list), numel(snr_list)); gn = zeros(numel(p_list), 1);
for ip = 1:numel(p_list)
  p = p_list(ip);
  for is = 0:numel(snr_list)                             % is = 0: noise only
    if is == 0, Aamp = 0; else, Aamp = 10^(snr_list(is)/20); end
    acc = 0; pn = 0; pn0 = 0;
    for t = 1:(T + (is == 0)*(Tn - T))
      nz = (randn(N,1) + 1j*randn(N,1))/sqrt(2); x = Aamp + nz; y = zeros(N,1); cnt = zeros(N,1);
      for dirn = [1 -1]
        if dirn > 0, tgt = (p+1:N)'; else, tgt = (1:N-p)'; end
        A = zeros(numel(tgt), p); for k = 1:p, A(:,k) = x(tgt - dirn*k); end
        G = A'*A; a = (G + mu*real(trace(G))/p*eye(p))\(A'*x(tgt)); y(tgt) = y(tgt) + A*a; cnt(tgt) = cnt(tgt) + 1;
      end
      y = y./cnt;
      if is == 0, pn = pn + abs(sum(y))^2; pn0 = pn0 + abs(sum(nz))^2; else, acc = acc + sum(y); end
    end
    if is == 0, gn(ip) = 10*log10(pn/pn0); else, gs(ip, is) = 20*log10(abs(acc/T)/(N*Aamp)); end
  end
end
fprintf('\nnoise-only DAS-output power gain (dB), %d trials, and forward-only estimate 10log10(p/(N-p)):\n', Tn);
for ip = 1:numel(p_list), fprintf('  p = %d: %6.1f dB   (estimate %6.1f dB)\n', p_list(ip), gn(ip), 10*log10(p_list(ip)/(N-p_list(ip)))); end
fprintf('echo amplitude gain (dB) at channel SNR %s dB:\n', sprintf('%4d', snr_list));
for ip = 1:numel(p_list), fprintf('  p = %d: %s\n', p_list(ip), sprintf('%6.1f', gs(ip,:))); end

%% 3. Image: plane-wave phantom (9 angles), DAS vs FXPF at -5 dB channel noise
fs = 20e6; xe = ((0:N-1) - (N-1)/2)*pitch; bw = 0.6; sigf = bw*f0/(2*sqrt(2*log(2))); sigt = 1/(2*pi*sigf);
ang = (-4:4)*pi/180; Na = numel(ang); t0 = 25e-6; Nt = 580; tt = t0 + (0:Nt-1)'/fs; noise_dB = -5; p = 4; W = 16;
nsc = round(13.6*8.6*80); xs = -6.8e-3 + 13.6e-3*rand(nsc,1); zs = 26e-3 + 8.6e-3*rand(nsc,1);
keep = (xs.^2 + (zs-30e-3).^2) > (2.5e-3)^2; xs = [xs(keep); -3e-3]; zs = [zs(keep); 24e-3]; as = [randn(sum(keep),1); 15];
IQ = zeros(Nt, N, Na); K4 = ceil(4*sigt*fs); E = repmat(1:N, numel(xs), 1);
for a = 1:Na
  tau = (xs*sin(ang(a)) + zs*cos(ang(a)))/c + sqrt((xs - xe).^2 + zs.^2)/c; ph = as.*exp(-1j*2*pi*f0*tau); n0 = (tau - t0)*fs;
  accr = zeros(Nt*N,1); acci = accr;
  for k = -K4:K4
    n = floor(n0) + k; g = exp(-((n - n0)/fs).^2/(2*sigt^2)); ok = n >= 0 & n < Nt; idx = n(ok) + 1 + (E(ok)-1)*Nt; v = g(ok).*ph(ok);
    accr = accr + accumarray(idx, real(v), [Nt*N 1]); acci = acci + accumarray(idx, imag(v), [Nt*N 1]);
  end
  IQ(:,:,a) = reshape(accr + 1j*acci, Nt, N);
end
Pref = mean(mean(abs(IQ(tt > 2*27.5e-3/c & tt < 2*33e-3/c, :, 5)).^2)); taps = exp(-((-K4:K4)'/fs).^2/(2*sigt^2));
for a = 1:Na, IQ(:,:,a) = IQ(:,:,a) + sqrt(Pref*10^(noise_dB/10)/sum(taps.^2))*conv2((randn(Nt+2*K4,N) + 1j*randn(Nt+2*K4,N))/sqrt(2), taps, 'valid'); end
dz = 0.075e-3; xg = (-6:0.15:6)*1e-3; zo = (22e-3:dz:34e-3)'; nzo = numel(zo); nx = numel(xg);
[Zp, Xp] = ndgrid(zo, xg); Xp = Xp(:); Zp = Zp(:); Ttx = (Xp*sin(ang) + Zp*cos(ang))/c; Ptx = exp(1j*2*pi*f0*Ttx); S = zeros(numel(Xp), N);
for e = 1:N
  trx = sqrt((Xp - xe(e)).^2 + Zp.^2)/c; acc = 0;
  for a = 1:Na
    fn = (Ttx(:,a) + trx - t0)*fs; n = floor(fn); fr = fn - n; ok = n >= 0 & n < Nt-1; col = IQ(:,e,a); v = zeros(size(fn));
    v(ok) = col(n(ok)+1).*(1-fr(ok)) + col(n(ok)+2).*fr(ok); acc = acc + v.*Ptx(:,a);
  end
  S(:,e) = acc.*exp(1j*2*pi*f0*trx);
end
S = reshape(S, nzo, nx, N);
% FXPF: per window, FFT along depth, per bin vectorized forward/backward LS over all columns, IFFT, overlap-add
tic; win = 0.5 - 0.5*cos(2*pi*(0:W-1)'/W); hop = W/2; Y = zeros(nzo, nx, N);   % periodic Hann: 50 % overlap sums to 1
for s0 = -hop:hop:nzo-1
  rr = s0 + (1:W); ok = rr >= 1 & rr <= nzo; seg = zeros(W, nx, N); seg(ok,:,:) = S(rr(ok),:,:);
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
      for k = 1:p, for i = k+1:p, fct = Gm(:,i,k)./Gm(:,k,k); for j = k:p, Gm(:,i,j) = Gm(:,i,j) - fct.*Gm(:,k,j); end, h(:,i) = h(:,i) - fct.*h(:,k); end, end
      cf_ = zeros(nx, p); for k = p:-1:1, sres = h(:,k); for j = k+1:p, sres = sres - Gm(:,k,j).*cf_(:,j); end, cf_(:,k) = sres./Gm(:,k,k); end
      P_ = zeros(nx, numel(tgt)); for q = 1:p, P_ = P_ + cf_(:,q).*A(:,:,q); end
      pr(:, tgt) = pr(:, tgt) + P_; cnt(tgt) = cnt(tgt) + 1;
    end
    G(kb,:,:) = reshape(pr./cnt, [1 nx N]);
  end
  g = ifft(G, [], 1); Y(rr(ok),:,:) = Y(rr(ok),:,:) + g(ok,:,:);
  if s0 == 72, Fchk = F; Gchk = G; end
end
tfx = toc;
imgs = {abs(sum(S,3))/N, abs(sum(Y,3))/N}; names = {'DAS', sprintf('FXPF p=%d W=%d', p, W)};
[ZI, XI] = ndgrid(zo, xg); cyst = (XI.^2 + (ZI-30e-3).^2) < (2e-3)^2; bgr = XI >= -5.8e-3 & XI <= -3.4e-3 & ZI >= 28e-3 & ZI <= 32e-3;
[~, rpt] = min(abs(zo - 24e-3)); sel = abs(xg + 3e-3) < 1.5e-3;
gcnr = @(a, b) 1 - sum(min(histc(a, linspace(min([a;b]), max([a;b]), 101))/numel(a), histc(b, linspace(min([a;b]), max([a;b]), 101))/numel(b)));
fprintf('\nimage at %d dB channel noise (FXPF %.1f s for %d x %d pixels)\n%-14s %8s %7s %7s %7s\n', noise_dB, tfx, nzo, nx, 'image', 'CR(dB)', 'gCNR', 'sSNR', 'pt -6dB');
for i = 1:2
  A = imgs{i}; pr = A(rpt, sel);
  fprintf('%-14s %8.1f %7.3f %7.2f %5.2fmm\n', names{i}, 10*log10(mean(A(bgr).^2)/mean(A(cyst).^2)), gcnr(A(cyst), A(bgr)), mean(A(bgr))/std(A(bgr),1), sum(pr >= max(pr)/2)*0.15);
end
% cross-check: the vectorized solver equals the per-vector solver of section 2 on one aperture vector
x = squeeze(Fchk(3, 40, :)); y = zeros(N,1); cnt = zeros(N,1);
for dirn = [1 -1]
  if dirn > 0, tgt = (p+1:N)'; else, tgt = (1:N-p)'; end
  A = zeros(numel(tgt), p); for k = 1:p, A(:,k) = x(tgt - dirn*k); end
  Gq = A'*A; a = (Gq + mu*real(trace(Gq))/p*eye(p))\(A'*x(tgt)); y(tgt) = y(tgt) + A*a; cnt(tgt) = cnt(tgt) + 1;
end
fprintf('vectorized vs per-vector prediction (window at row 73, bin 3, column 40): max rel. difference %.1e\n', max(abs(y./cnt - squeeze(Gchk(3,40,:))))/max(abs(y)));

%% Figures
figure('Position', [50 50 1150 330]);
subplot(1,3,1); hold on; cols = 'brmk';
for ip = 1:numel(p_list), plot(snr_list, gs(ip,:), [cols(ip) '-o']); plot(snr_list([1 end]), gn(ip)*[1 1], [cols(ip) ':']); end
xlabel('channel SNR (dB)'); ylabel('DAS-output gain (dB)'); title('FXPF: echo (solid) and noise (dotted)'); legend('p=1','','p=2','','p=4','','p=8',''); grid on;
for i = 1:2
  subplot(1,3,i+1); A = imgs{i}/max(imgs{i}(:)); imagesc(xg*1e3, zo*1e3, 20*log10(A + 1e-12), [-60 0]); axis image; colormap(gray);
  title([names{i} ', ' num2str(noise_dB) ' dB noise']); xlabel('x (mm)'); ylabel('z (mm)');
end
