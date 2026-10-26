%% m25_dmas.m  -  delay-multiply-and-sum (DMAS / F-DMAS)
% Purpose : verify the O(N) identity sum_{i<j} s_i s_j = ((sum s)^2 - sum s^2)/2 for the
%           signed-square-root products, show the DC + 2f0 spectrum of the DMAS sum,
%           compute the narrowband sum-coarray beam (DMAS at 2f0) against DAS, and
%           compare DAS, unfiltered DMAS and F-DMAS images (points, cyst in speckle).
% Change  : noise_dB, the band-pass width (sigh), N.
% Expect  : identity error ~1e-13; narrowband PSL -13.3 dB (DAS) vs -26.5 dB (DMAS),
%           width ratio 0.74; images: F-DMAS PSL ~-27 dB vs ~-19 dB DAS.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8. Runtime ~6 s.
clear; close all; rng(4);
c = 1540; f0 = 5e6; fs = 20e6; lam = c/f0; N = 64; pitch = 0.3e-3; xe = ((0:N-1) - (N-1)/2)*pitch;
bw = 0.6; sigf = bw*f0/(2*sqrt(2*log(2))); sigt = 1/(2*pi*sigf); ang = (-4:4)*pi/180; Na = numel(ang);

%% 1. The O(N) identity (random RF aperture vectors)
s = randn(1000, N); sh = sign(s).*sqrt(abs(s)); slow = zeros(1000,1);
for i = 1:N-1, for j = i+1:N, slow = slow + sign(s(:,i).*s(:,j)).*sqrt(abs(s(:,i).*s(:,j))); end, end
fast = (sum(sh,2).^2 - sum(abs(s),2))/2;
fprintf('identity: max |pairwise - O(N)| = %.2e (relative %.1e); %d pair products replaced by %d square roots\n', max(abs(slow - fast)), max(abs(slow - fast))/max(abs(slow)), N*(N-1)/2, N);

%% 2. Spectrum of the DMAS sum for an aligned echo: every pair term equals |rf(t)|
fsr = 200e6; t = (-2048:2047)'/fsr; rf = exp(-t.^2/(2*sigt^2)).*cos(2*pi*f0*t); Y = abs(fft(abs(rf))); fsp = (0:4095)'*fsr/4096;
e0 = sum(Y(fsp < 2.5e6).^2); e2 = sum(Y(fsp > 7.5e6 & fsp < 12.5e6).^2); e4 = sum(Y(fsp > 17.5e6 & fsp < 22.5e6).^2);
fprintf('aligned echo, |rf| spectrum: energy DC band : 2f0 band : 4f0 band = 1 : %.3f : %.4f (positive frequencies)\n', e2/e0, e4/e0);

%% 3. Narrowband sum coarray at z = 30 mm: 2f0 term of pair (i,j) has phase phi_i + phi_j
x = linspace(-8e-3, 8e-3, 32001); q = 2*pi*pitch*x/(lam*30e-3); AF = exp(1j*q'*(0:N-1)); a1 = sum(AF,2); a2 = sum(AF.^2,2);
dasNB = abs(a1)/N; dmNB = abs((a1.^2 - a2)/2)/(N*(N-1)/2);
wid = @(y) max(x(y >= 0.5)) - min(x(y >= 0.5));
i0 = 16001; il = i0; while dasNB(il-1) <= dasNB(il), il = il - 1; end
jl = i0; while dmNB(jl-1) <= dmNB(jl), jl = jl - 1; end
fprintf('narrowband: DAS FWHM %.3f mm PSL %.1f dB; DMAS(2f0) FWHM %.3f mm PSL %.1f dB; width ratio %.3f\n', wid(dasNB)*1e3, 20*log10(max(dasNB(1:il))), wid(dmNB)*1e3, 20*log10(max(dmNB(1:jl))), wid(dmNB)/wid(dasNB));

%% 4. Plane-wave simulation: points at z = 24 mm, cyst in speckle, noise
t0 = 25e-6; Nt = 580; tt = t0 + (0:Nt-1)'/fs; noise_dB = -20;
nsc = round(13.6*8.6*80); xs = (rand(nsc,1)*13.6 - 6.8)*1e-3; zs = (26 + 8.6*rand(nsc,1))*1e-3;
keep = (xs.^2 + (zs-30e-3).^2) > (2.5e-3)^2; xs = [xs(keep); -3e-3; 2.25e-3; 2.75e-3]; zs = [zs(keep); 24e-3; 24e-3; 24e-3]; as = [randn(sum(keep),1); 15; 15; 15];
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

%% 5. RF aperture data on a fine axial grid (lambda/16), then DAS, DMAS, F-DMAS
sub = 4; dz = 0.075e-3; dzf = dz/sub; fsz = c/(2*dzf); xg = (-6:0.05:6)*1e-3; zf = (22e-3:dzf:34e-3)';
[Zp, Xp] = ndgrid(zf, xg); Xp = Xp(:); Zp = Zp(:); S = zeros(numel(Xp), N); Ttx = (Xp*sin(ang) + Zp*cos(ang))/c; Ptx = exp(1j*2*pi*f0*Ttx);
for e = 1:N
  trx = sqrt((Xp - xe(e)).^2 + Zp.^2)/c; acc = 0;
  for a = 1:Na
    fn = (Ttx(:,a) + trx - t0)*fs; n = floor(fn); fr = fn - n; ok = n >= 0 & n < Nt-1; col = IQ(:,e,a); v = zeros(size(fn));
    v(ok) = col(n(ok)+1).*(1-fr(ok)) + col(n(ok)+2).*fr(ok); acc = acc + v.*Ptx(:,a);
  end
  S(:,e) = acc.*exp(1j*2*pi*f0*trx);
end
nzf = numel(zf); nx = numel(xg); R = reshape(real(S), nzf, nx, N); Sc = reshape(S, nzf, nx, N);
tic; das = abs(sum(Sc,3))/N; tdas = toc;
tic; shh = sign(R).*sqrt(abs(R)); yd = (sum(shh,3).^2 - sum(abs(R),3))/2;
nf = 2^nextpow2(nzf + 64); f = [0:nf/2, -nf/2+1:-1]'*fsz/nf; sigh = sqrt(2)*sigf;
Yd = fft(yd, nf, 1); envF = abs(ifft(Yd.*(2*exp(-(f - 2*f0).^2/(2*sigh^2)).*(f > 0)), [], 1)); envF = envF(1:nzf,:);
envU = abs(ifft(Yd.*(2*(f > 0) + (f == 0)), [], 1)); envU = envU(1:nzf,:); tdm = toc;   % unfiltered: analytic-signal envelope
imgs = {das(1:sub:end,:), envU(1:sub:end,:), envF(1:sub:end,:)}; names = {'DAS', 'DMAS, no filter', 'F-DMAS'}; zo = zf(1:sub:end);
fprintf('time per image: DAS %.2f s, DMAS + band-pass %.2f s (%d x %d RF samples)\n', tdas, tdm, nzf, nx);
[ZI, XI] = ndgrid(zo, xg); cyst = (XI.^2 + (ZI-30e-3).^2) < (2e-3)^2; bgr = XI >= -5.8e-3 & XI <= -3.4e-3 & ZI >= 28e-3 & ZI <= 32e-3;
[~, rp] = min(abs(zo - 24e-3)); xsel = xg >= -5.5e-3 & xg <= -0.5e-3; xq = xg(xsel);
fprintf('\n%-16s %6s %6s %6s %6s %6s\n', 'image', 'FWHM', 'PSL', 'CR', 'gCNR', 'sSNR');
for i = 1:3
  A = imgs{i}; pr = A(rp, xsel); [pk, ip] = max(pr); half = pk/2;
  il = ip; while pr(il) > half, il = il - 1; end; ir = ip; while pr(ir) > half, ir = ir + 1; end
  xl = xq(il) + (half - pr(il))/(pr(il+1) - pr(il))*(xq(il+1) - xq(il)); xr = xq(ir-1) + (half - pr(ir-1))/(pr(ir) - pr(ir-1))*(xq(ir) - xq(ir-1));
  ml = il; while ml > 1 && pr(ml-1) <= pr(ml), ml = ml - 1; end; mr = ir; while mr < numel(pr) && pr(mr+1) <= pr(mr), mr = mr + 1; end
  side = pr; side(ml:mr) = 0; ci = A(cyst); bi = A(bgr); v = [ci; bi]; lab = [zeros(numel(ci),1); ones(numel(bi),1)]; [~, o] = sort(v); lab = lab(o);
  bins = floor((0:numel(v)-1)'*100/numel(v)); g = 1 - sum(min(accumarray(bins(lab==0)+1,1,[100 1])/numel(ci), accumarray(bins(lab==1)+1,1,[100 1])/numel(bi)));
  fprintf('%-16s %5.3fmm %6.1f %6.1f %6.3f %6.2f\n', names{i}, (xr - xl)*1e3, 20*log10(max(side)/pk), 10*log10(mean(bi.^2)/mean(ci.^2)), g, mean(bi)/std(bi,1));
end
fprintf('(FWHM from the 0.05 mm image columns through the point at x = -3 mm, z = 24 mm)\n');

%% Figures
figure('Position', [50 50 1200 330]);
for i = 1:3
  subplot(1,4,i); A = imgs{i}/max(imgs{i}(:)); imagesc(xg*1e3, zo*1e3, 20*log10(A + 1e-12), [-60 0]); axis image; colormap(gray); title(names{i}); xlabel('x (mm)'); ylabel('z (mm)');
end
subplot(1,4,4); plot(fsp(1:500)/1e6, Y(1:500)/Y(1), 'b'); hold on; plot(fsp(1:500)/1e6, 0.5*exp(-(fsp(1:500) - 2*f0).^2/(2*sigh^2)), 'r--');
xlabel('frequency (MHz)'); ylabel('|spectrum| (norm.)'); title('DMAS sum of an aligned echo'); legend('|rf| spectrum', 'band-pass');
figure; plot(x*1e3, 20*log10(dasNB + 1e-9), 'k', x*1e3, 20*log10(dmNB + 1e-9), 'b'); ylim([-60 0]); xlim([-4 4]);
xlabel('lateral position (mm)'); ylabel('beam (dB)'); legend('DAS', 'DMAS at 2f_0 (sum coarray)'); title('Narrowband: DMAS = triangular aperture at \lambda/2');
