%% m25_slsc.m  -  short-lag spatial coherence (SLSC) imaging
% Purpose : simulate 64-channel plane-wave data (9 angles, -4..4 deg) of speckle with an
%           anechoic cyst and a point target, measure the spatial coherence R(m) of the
%           delayed RF in speckle, at the point and in the cyst, compare speckle R(m)
%           with the van Cittert-Zernike (VCZ) triangle of the synthesized transmit
%           aperture, show the noise scaling SNR/(1+SNR), and sweep the maximum lag M.
% Change  : noise_dB (channel noise re speckle), M_list, Ks (axial kernel, fine samples).
% Expect  : speckle R(m) close to 1 - m/m_tx with m_tx = 15.7 lags; R_sl(M=8) ~ 5.7 (VCZ);
%           at 0 dB channel SNR the lags m >= 1 halve; point R(m) ~ 1; cyst R(m >= 1) ~ 0;
%           larger M -> narrower point, lower speckle SNR. Runtime ~15 s in Octave.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(2);

%% Acquisition, phantom, channel data (signal and unit noise kept separate: beamforming is linear)
c = 1540; f0 = 5e6; fs = 20e6; lam = c/f0; N = 64; pitch = 0.3e-3;        % m/s, Hz, Hz, m, -, m
xe = ((0:N-1) - (N-1)/2)*pitch; bw = 0.6; sigf = bw*f0/(2*sqrt(2*log(2))); sigt = 1/(2*pi*sigf);
ang = (-4:4)*pi/180; Na = numel(ang); t0 = 25e-6; Nt = 580; tt = t0 + (0:Nt-1)'/fs;
noise_dB = -20; Ks = 8; M_list = [1 2 4 8 12 16 24];
nsc = round(13.6*8.6*80); xs = -6.8e-3 + 13.6e-3*rand(nsc,1); zs = 26e-3 + 8.6e-3*rand(nsc,1);
keep = (xs.^2 + (zs-30e-3).^2) > (2.5e-3)^2; xs = [xs(keep); -3e-3]; zs = [zs(keep); 24e-3]; as = [randn(sum(keep),1); 15];
tic; IQ = zeros(Nt, N, Na); K4 = ceil(4*sigt*fs); E = repmat(1:N, numel(xs), 1);
for a = 1:Na
  tau = (xs*sin(ang(a)) + zs*cos(ang(a)))/c + sqrt((xs - xe).^2 + zs.^2)/c; ph = as.*exp(-1j*2*pi*f0*tau); n0 = (tau - t0)*fs;
  accr = zeros(Nt*N,1); acci = accr;
  for k = -K4:K4
    n = floor(n0) + k; g = exp(-((n - n0)/fs).^2/(2*sigt^2)); ok = n >= 0 & n < Nt; idx = n(ok) + 1 + (E(ok)-1)*Nt; v = g(ok).*ph(ok);
    accr = accr + accumarray(idx, real(v), [Nt*N 1]); acci = acci + accumarray(idx, imag(v), [Nt*N 1]);
  end
  IQ(:,:,a) = reshape(accr + 1j*acci, Nt, N);
end
Pref = mean(mean(abs(IQ(tt > 2*27.5e-3/c & tt < 2*33e-3/c, :, 5)).^2));   % mean speckle channel power
taps = exp(-((-K4:K4)'/fs).^2/(2*sigt^2)); NQ = zeros(Nt, N, Na);
for a = 1:Na, NQ(:,:,a) = sqrt(Pref/sum(taps.^2))*conv2((randn(Nt+2*K4,N) + 1j*randn(Nt+2*K4,N))/sqrt(2), taps, 'valid'); end
fprintf('simulated %d scatterers in %.1f s\n', numel(xs), toc);

%% Delayed, compounded aperture data on a fine axial grid (dz/4, RF rate 41 MHz)
sub = 4; dz = 0.075e-3; dzf = dz/sub; xg = (-6:0.15:6)*1e-3; zf = (22e-3:dzf:34e-3)'; nzf = numel(zf); nx = numel(xg);
[Zp, Xp] = ndgrid(zf, xg); Xp = Xp(:); Zp = Zp(:); Ttx = (Xp*sin(ang) + Zp*cos(ang))/c; Ptx = exp(1j*2*pi*f0*Ttx);
Ssig = zeros(numel(Xp), N); Snoi = Ssig;
for e = 1:N
  trx = sqrt((Xp - xe(e)).^2 + Zp.^2)/c; as_ = 0; an_ = 0;
  for a = 1:Na
    fn = (Ttx(:,a) + trx - t0)*fs; n = floor(fn); fr = fn - n; ok = n >= 0 & n < Nt-1;
    v = zeros(size(fn)); col = IQ(:,e,a); v(ok) = col(n(ok)+1).*(1-fr(ok)) + col(n(ok)+2).*fr(ok); as_ = as_ + v.*Ptx(:,a);
    v = zeros(size(fn)); col = NQ(:,e,a); v(ok) = col(n(ok)+1).*(1-fr(ok)) + col(n(ok)+2).*fr(ok); an_ = an_ + v.*Ptx(:,a);
  end
  Ssig(:,e) = as_.*exp(1j*2*pi*f0*trx); Snoi(:,e) = an_.*exp(1j*2*pi*f0*trx);
end
Ssig = reshape(Ssig, nzf, nx, N); Snoi = reshape(Snoi, nzf, nx, N);
ro = (1:sub:nzf)'; nzo = numel(ro); zo = zf(ro); lo = max(1, ro - Ks); hi = min(nzf, ro + Ks);
[ZI, XI] = ndgrid(zo, xg); cyst = (XI.^2 + (ZI-30e-3).^2) < (2e-3)^2; bgr = XI >= -5.8e-3 & XI <= -3.4e-3 & ZI >= 28e-3 & ZI <= 32e-3;

%% 1. Coherence curves R(m), m = 0..40, in speckle (noise-free and 0 dB), at the point, in the cyst
m_tx = Na*30e-3*(sin(ang(2)) - sin(ang(1)))/pitch;                        % VCZ: synthesized transmit width in lags
fprintf('VCZ: 9 angles x 1 deg at z = 30 mm synthesize %.2f mm = %.1f element lags\n', m_tx*pitch*1e3, m_tx);
cases = {Inf, 0}; Rc = zeros(41, 4); ipt = find(abs(xg + 3e-3) < 1e-9); [~, rpt] = min(abs(zo - 24e-3));
for ic = 1:2
  snr = cases{ic}; if snr == Inf, R = real(Ssig); else, R = real(Ssig + 10^(-snr/20)*Snoi); end
  Cs = cumsum(cat(1, zeros(1, nx, N), R.^2), 1); En = Cs(hi+1,:,:) - Cs(lo,:,:);
  for m = 0:40
    Cs = cumsum(cat(1, zeros(1, nx, N-m), R(:,:,1:N-m).*R(:,:,1+m:N)), 1); num = Cs(hi+1,:,:) - Cs(lo,:,:);
    Rm = mean(num./sqrt(En(:,:,1:N-m).*En(:,:,1+m:N) + eps), 3);
    Rc(m+1, ic) = mean(Rm(bgr));
    if ic == 1, Rc(m+1, 3) = Rm(rpt, ipt); end
    if ic == 2, Rc(m+1, 4) = mean(Rm(cyst)); end
  end
end
vcz = max(0, 1 - (0:40)'/m_tx);
fprintf('speckle R(1), R(4), R(8), R(16): sim %.2f %.2f %.2f %.2f | VCZ %.2f %.2f %.2f %.2f\n', Rc([2 5 9 17],1), vcz([2 5 9 17]));
fprintf('R_sl(M = 8) in speckle: sim %.2f, VCZ %.2f; at 0 dB channel SNR: sim %.2f, predicted x SNR/(1+SNR) = %.2f\n', ...
  sum(Rc(2:9,1)), sum(vcz(2:9)), sum(Rc(2:9,2)), 0.5*sum(Rc(2:9,1)));
fprintf('point (-3, 24) mm: R(1) = %.3f, R_sl(8) = %.2f (ideal 8);  cyst at 0 dB SNR: mean R(1..8) = %.3f\n', Rc(2,3), sum(Rc(2:9,3)), mean(Rc(2:9,4)));

%% 2. Images: DAS and SLSC at the default noise level, sweep of the maximum lag M
R = real(Ssig + 10^(noise_dB/20)*Snoi); Sc = Ssig(ro,:,:) + 10^(noise_dB/20)*Snoi(ro,:,:);
das = abs(sum(Sc,3))/N; Cs = cumsum(cat(1, zeros(1, nx, N), R.^2), 1); En = Cs(hi+1,:,:) - Cs(lo,:,:);
Rlag = zeros(nzo, nx, max(M_list));
for m = 1:max(M_list)
  Cs = cumsum(cat(1, zeros(1, nx, N-m), R(:,:,1:N-m).*R(:,:,1+m:N)), 1); num = Cs(hi+1,:,:) - Cs(lo,:,:);
  Rlag(:,:,m) = mean(num./sqrt(En(:,:,1:N-m).*En(:,:,1+m:N) + eps), 3);
end
gcnr = @(a, b) 1 - sum(min(histc(a, linspace(min([a;b]), max([a;b]), 101))/numel(a), histc(b, linspace(min([a;b]), max([a;b]), 101))/numel(b)));
fw = @(pr) sum(pr >= max(pr)/2)*0.15;                                     % -6 dB width on the 0.15 mm grid (mm)
met = zeros(numel(M_list)+1, 4);
A = das; met(1,:) = [fw(A(rpt, abs(xg + 3e-3) < 2.6e-3)), 10*log10(mean(A(bgr).^2)/mean(A(cyst).^2)), gcnr(A(cyst), A(bgr)), mean(A(bgr))/std(A(bgr),1)];
for k = 1:numel(M_list)
  A = max(sum(Rlag(:,:,1:M_list(k)), 3), 0); ci = max(A(cyst), 1e-6); bi = A(bgr);
  met(k+1,:) = [fw(A(rpt, abs(xg + 3e-3) < 2.6e-3)), 10*log10(mean(bi.^2)/mean(ci.^2)), gcnr(ci, bi), mean(bi)/std(bi,1)];
end
fprintf('\n%-10s %10s %8s %7s %8s\n', 'image', 'point -6dB', 'CR(dB)', 'gCNR', 'sSNR');
fprintf('%-10s %8.2fmm %8.1f %7.3f %8.2f\n', 'DAS', met(1,:));
for k = 1:numel(M_list), fprintf('SLSC M=%-3d %8.2fmm %8.1f %7.3f %8.2f\n', M_list(k), met(k+1,:)); end
fprintf('(point width is quantized to the 0.15 mm pixel grid; SLSC enters the metrics as max(R_sl, 0))\n');

%% Figures
figure('Position', [50 50 1200 380]);
subplot(1,3,1); imagesc(xg*1e3, zo*1e3, 20*log10(das/max(das(:)) + 1e-12), [-60 0]); axis image; colormap(gray); title('DAS (60 dB)'); xlabel('x (mm)'); ylabel('z (mm)');
A = max(sum(Rlag(:,:,1:8), 3), 0); subplot(1,3,2); imagesc(xg*1e3, zo*1e3, A/max(A(:)), [0 1]); axis image; title('SLSC, M = 8 (linear)'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(1,3,3); plot(0:40, Rc(:,1), 'b', 0:40, Rc(:,2), 'b--', 0:40, vcz, 'k:', 0:40, Rc(:,3), 'r', 0:40, Rc(:,4), 'm', 'LineWidth', 1.3);
xlabel('element lag m'); ylabel('R(m)'); title('Spatial coherence'); legend('speckle', 'speckle, 0 dB SNR', 'VCZ triangle', 'point', 'cyst (noise)'); ylim([-0.2 1.05]); grid on;
figure('Position', [50 50 1100 300]);
lab = {'point -6 dB width (mm)', 'CR (dB)', 'gCNR', 'speckle SNR'};
for q = 1:4, subplot(1,4,q); plot(M_list, met(2:end,q), 'ro-', M_list, met(1,q)*ones(size(M_list)), 'k--'); xlabel('maximum lag M'); title(lab{q}); grid on; end
