%% m25_coherence_factors.m  -  CF, Wiener postfilter, GCF, PCF and SCF on one dataset
% Purpose : compute the five coherence weights from the same plane-wave aperture data
%           (cyst in speckle, bright off-axis scatterer, channel noise), verify the
%           identity H = N CF/((N-1) CF + 1), compare the speckle CF with the VCZ
%           prediction for this transmit sequence, and report CR, gCNR, speckle SNR.
% Change  : noise_dB, M0 list, gam (PCF), pSCF.
% Expect  : speckle CF ratio of means ~0.23 (VCZ triangle 0.225), CF in noise 1/N;
%           CR rises CF > PCF > GCF > Wiener > DAS while gCNR does not follow.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8. Runtime ~6 s.
clear; close all; rng(3);
c = 1540; f0 = 5e6; fs = 20e6; lam = c/f0; N = 64; pitch = 0.3e-3; xe = ((0:N-1) - (N-1)/2)*pitch;
bw = 0.6; sigf = bw*f0/(2*sqrt(2*log(2))); sigt = 1/(2*pi*sigf); ang = (-4:4)*pi/180; Na = numel(ang);
t0 = 25e-6; Nt = 580; tt = t0 + (0:Nt-1)'/fs; noise_dB = -15;

%% Simulation: speckle with an anechoic cyst (r = 2.5 mm) and a bright scatterer 4.2 mm off-axis
nsc = round(13.6*8.6*80); xs = (rand(nsc,1)*13.6 - 6.8)*1e-3; zs = (26 + 8.6*rand(nsc,1))*1e-3;
keep = (xs.^2 + (zs-30e-3).^2) > (2.5e-3)^2; xs = [xs(keep); 4.2e-3]; zs = [zs(keep); 30.5e-3]; as = [randn(sum(keep),1); 60];
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

%% Aperture data on the image grid
xg = (-6:0.15:6)*1e-3; zg = (26:0.075:34)*1e-3; [Zp, Xp] = ndgrid(zg, xg); Xp = Xp(:); Zp = Zp(:); S = zeros(numel(Xp), N);
Ttx = (Xp*sin(ang) + Zp*cos(ang))/c; Ptx = exp(1j*2*pi*f0*Ttx);
for e = 1:N
  trx = sqrt((Xp - xe(e)).^2 + Zp.^2)/c; acc = 0;
  for a = 1:Na
    fn = (Ttx(:,a) + trx - t0)*fs; n = floor(fn); fr = fn - n; ok = n >= 0 & n < Nt-1; col = IQ(:,e,a); v = zeros(size(fn));
    v(ok) = col(n(ok)+1).*(1-fr(ok)) + col(n(ok)+2).*fr(ok); acc = acc + v.*Ptx(:,a);
  end
  S(:,e) = acc.*exp(1j*2*pi*f0*trx);
end
nz = numel(zg); nx = numel(xg);

%% Weights
ybar = mean(S, 2); P = mean(abs(S).^2, 2); sig2 = mean(abs(S - ybar).^2, 2);
CF = abs(sum(S,2)).^2./(N*sum(abs(S).^2,2) + eps);                  % |sum s|^2/(N sum |s|^2)
Hw = abs(ybar).^2./(abs(ybar).^2 + sig2/N + eps);                    % Wiener postfilter from its definition
fprintf('identity check: max |H - N CF/((N-1)CF+1)| = %.2e,  max |P - (|ybar|^2 + sig2)| = %.2e\n', max(abs(Hw - N*CF./((N-1)*CF + 1))), max(abs(P - abs(ybar).^2 - sig2))/max(P));
Sk = fft(S, [], 2); M0s = [0 1 3]; GCF = zeros(numel(Xp), numel(M0s));
for i = 1:numel(M0s), low = abs(Sk(:,1)).^2; for k = 1:M0s(i), low = low + abs(Sk(:,1+k)).^2 + abs(Sk(:,N+1-k)).^2; end, GCF(:,i) = low./(sum(abs(Sk).^2,2) + eps); end
fprintf('GCF(M0 = 0) equals CF: max difference %.2e\n', max(abs(GCF(:,1) - CF)));
gam = 1; phs = angle(S); pha = phs - pi*sign(phs); PCF = max(0, 1 - gam*min(std(phs,1,2), std(pha,1,2))/(pi/sqrt(3)));
pSCF = 1; bbar = mean(sign(real(S)), 2); SCF = abs(1 - sqrt(max(0, 1 - bbar.^2))).^pSCF;

%% Statistics: speckle CF vs VCZ, noise CF vs 1/N
[ZI, XI] = ndgrid(zg, xg); cyst = (XI.^2 + (ZI-30e-3).^2) < (2e-3)^2; bgr = XI >= -5.8e-3 & XI <= -3.4e-3 & ZI >= 28e-3 & ZI <= 32e-3;
mtx = Na*30e-3*(sin(ang(2)) - sin(ang(1)))/pitch; m = -(N-1):(N-1); Ctri = max(0, 1 - abs(m)/mtx);
Sb = S(bgr(:), :);
fprintf('speckle: CF ratio of means %.3f (VCZ triangle, m_tx = %.1f: %.3f), mean of per-pixel CF %.3f\n', mean(abs(sum(Sb,2)).^2)/(N*mean(sum(abs(Sb).^2,2))), mtx, sum((N - abs(m)).*Ctri)/N^2, mean(CF(bgr(:))));
Wn = (randn(20000, N) + 1j*randn(20000, N))/sqrt(2); cfn = abs(sum(Wn,2)).^2./(N*sum(abs(Wn).^2,2));
fprintf('white noise: mean CF %.4f (1/N = %.4f), mean Wiener H %.3f (N/(2N-1) = %.3f at CF = 1/N), mean PCF %.3f\n', mean(cfn), 1/N, mean(N*cfn./((N-1)*cfn + 1)), N/(2*N-1), mean(max(0, 1 - std(angle(Wn),1,2)/(pi/sqrt(3)))));

%% Image metrics (amplitude = weight x |DAS|)
das = abs(ybar); imgs = [das, CF.*das, Hw.*das, GCF(:,2).*das, GCF(:,3).*das, PCF.*das, SCF.*das];
names = {'DAS', 'CF', 'Wiener', 'GCF M0=1', 'GCF M0=3', 'PCF', 'SCF'};
fprintf('\n%-9s %6s %6s %6s %8s\n', 'method', 'CR', 'gCNR', 'sSNR', 'w(spk)');
for i = 1:numel(names)
  A = imgs(:,i); ci = max(A(cyst(:)), 1e-12*max(A)); bi = A(bgr(:)); v = [ci; bi]; lab = [zeros(numel(ci),1); ones(numel(bi),1)];
  [vs, o] = sort(v); lab = lab(o); bins = floor((0:numel(v)-1)'*100/numel(v)); [~, first, ic] = unique(vs, 'first'); bins = bins(first(ic));
  g = 1 - sum(min(accumarray(bins(lab==0)+1,1,[100 1])/numel(ci), accumarray(bins(lab==1)+1,1,[100 1])/numel(bi)));
  fprintf('%-9s %6.1f %6.3f %6.2f %8.3f\n', names{i}, 10*log10(mean(bi.^2)/mean(ci.^2)), g, mean(bi)/std(bi,1), mean(A(bgr(:))./das(bgr(:))));
end

%% Figures
figure('Position', [50 50 1300 330]);
for i = 1:numel(names)
  subplot(1, numel(names), i); A = reshape(imgs(:,i), nz, nx); A = A/max(A(:));
  imagesc(xg*1e3, zg*1e3, 20*log10(A + 1e-12), [-60 0]); axis image; colormap(gray); title(names{i}); xlabel('x (mm)'); ylabel('z (mm)');
end
figure; edges = linspace(0, 1, 41); hs = histc(CF(bgr(:)), edges); hn = histc(cfn, edges);
bar(edges + 0.0125, [hs/sum(hs) hn/sum(hn)]*40, 1); xlim([0 1]); legend('speckle (background ROI)', 'white noise'); xlabel('CF'); ylabel('pdf'); title('CF is broadly distributed in speckle');
