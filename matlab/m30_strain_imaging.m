%% m30_strain_imaging.m
% Purpose : 1D quasi-static strain imaging. Pre- and post-compression RF
%           lines are synthesized from 2000 point scatterers (5 MHz Gaussian
%           pulse, 60 % bandwidth, fs = 40 MHz). A stiff layer (18-24 mm,
%           modulus ratio R) sits in a uniform-stress (series spring) medium.
%           Displacement by normalized cross-correlation with a guided search
%           and parabolic subsample fit; strain by least squares over N
%           estimates. Then the strain filter: SNRe versus applied strain in a
%           homogeneous medium, with and without companding (stretching the
%           post-compression line by the mean strain).
% Change  : epsApp, R, W (kernel, m), N (LSQ span), snrdB, the strain list.
% Expect  : at 1 % strain, 2 mm kernel, N = 9, 30 dB: background strain
%           1.111 %, layer 0.370 %, measured strain ratio close to 3, SNRe of
%           order 5 (speckle-limited); companding raises SNRe at large strain.
%           Eshelby cylinder: strain ratio (1+r)/2 = 2.5 and 5.5 for r = 4, 10.
rng(1);

%% Parameters
c = 1540; fs = 40e6; f0 = 5e6; B = 0.6; Z = 0.04; nsc = 2000;
dzs = c/(2*fs); ns = round(Z/dzs);               % 19.25 um per sample
sf = B*f0/2/sqrt(2*log(2)); st = 1/(2*pi*sf);    % pulse envelope std (s)
half = ceil(4*st*fs);
zsc = rand(nsc,1)*Z; amp = randn(nsc,1);
n1 = randn(ns,1); n2 = randn(ns,1);
z1 = 0.018; z2 = 0.024; w = z2 - z1;

% RF synthesis (vectorized with accumarray)

%% Run one configuration
W = 2e-3; N = 9; snrdB = 30; epsApp = 0.01; R = 3;
ebg = epsApp*Z/(Z - w + w/R); eL = ebg/R;
ufun = @(z) -( (z<z1).*ebg.*z + (z>=z1 & z<z2).*(ebg*z1 + eL*(z-z1)) + (z>=z2).*(ebg*z1 + eL*w + ebg*(z-z2)) );
fprintf('True strain: background %.3f %%, layer %.3f %% (ratio %.2f)\n', ebg*100, eL*100, ebg/eL);
fprintf('Max displacement %.0f um = %.1f samples\n', -ufun(Z)*1e6, -ufun(Z)/dzs);

idx0 = (-half:half);
synth = @(zk) accumarray(reshape(min(max(round(2*zk/c*fs) + idx0, 1), ns), [], 1), ...
          reshape(repmat(amp, 1, numel(idx0)) .* exp(-((min(max(round(2*zk/c*fs) + idx0, 1), ns) - 1)/fs - 2*zk/c).^2/(2*st^2)) ...
          .* cos(2*pi*f0*((min(max(round(2*zk/c*fs) + idx0, 1), ns) - 1)/fs - 2*zk/c)), [], 1), [ns 1]);
addn = @(y, nz, snr) y + sqrt(mean(y.^2)/10^(snr/10))*nz;

x = addn(synth(zsc), n1, snrdB);
y = addn(synth(zsc + ufun(zsc)), n2, snrdB);

nw = round(W/dzs); step = 13; S = 12;
zc = []; uu = []; prev = 0;
for c0 = 41:step:(ns - nw - 41)
  xs = x(c0:c0+nw-1); ex = sum(xs.^2); lags = (prev-S):(prev+S); cc = -2*ones(size(lags));
  for q = 1:numel(lags)
    i0 = c0 + lags(q); if i0 < 1 || i0 + nw - 1 > ns, continue; end
    ys = y(i0:i0+nw-1); cc(q) = sum(xs.*ys)/sqrt(ex*sum(ys.^2));
  end
  [~, q] = max(cc); d = 0;
  if q > 1 && q < numel(lags), den = cc(q-1) - 2*cc(q) + cc(q+1); if den < 0, d = 0.5*(cc(q-1) - cc(q+1))/den; end, end
  prev = lags(q); zc(end+1) = (c0 - 1 + nw/2)*dzs; uu(end+1) = (lags(q) + d)*dzs; %#ok<AGROW>
end
ne = numel(zc) - N + 1; ze = zeros(1,ne); ee = zeros(1,ne);
for i = 1:ne
  zz = zc(i:i+N-1); u2 = uu(i:i+N-1); zm = mean(zz);
  ee(i) = -sum((zz - zm).*(u2 - mean(u2)))/sum((zz - zm).^2); ze(i) = zm;
end
bg = ze > 0.004 & ze < 0.016; bg2 = ze > 0.028 & ze < 0.036; ly = ze > 0.0195 & ze < 0.0225;
SNRe = mean(ee(bg))/std(ee(bg));
fprintf('Kernel %.1f mm, N = %d, %d dB: SNRe (4-16 mm) = %.1f, strain std = %.3f %%\n', W*1e3, N, snrdB, SNRe, std(ee(bg))*100);
fprintf('Measured strain ratio = %.2f (true %.1f)\n', mean([ee(bg) ee(bg2)])/mean(ee(ly)), R);
fprintf('Eshelby cylinder strain ratio (1+r)/2: r = 4 -> %.1f, r = 10 -> %.1f\n', (1+4)/2, (1+10)/2);

figure(1);
subplot(2,1,1); zt = linspace(0, Z, 400);
plot(zt*1e3, ufun(zt)*1e6, 'k--', zc*1e3, uu*1e6, 'b'); ylabel('u (\mum)'); title('Displacement: true (dashed) and NCC estimate');
subplot(2,1,2);
plot([0 18 18 24 24 40], [ebg ebg eL eL ebg ebg]*100, 'k--', ze*1e3, ee*100, 'b'); xlabel('depth (mm)'); ylabel('strain (%)'); title('LSQ strain');

%% Strain filter (homogeneous medium)
epsList = [0.0005 0.001 0.003 0.01 0.03 0.06 0.1];
SN = zeros(2, numel(epsList));
for ie = 1:numel(epsList)
  e = epsList(ie);
  yh = addn(synth(zsc*(1 - e)), n2, snrdB);
  for comp = 0:1
    yy = yh; sh = 0;
    if comp, sh = e; yy = interp1((0:ns-1)', yh, (0:ns-1)'*(1 - sh), 'linear', 0); end
    zc = []; uu = []; prev = 0;
    for c0 = 41:step:(ns - nw - 41)
      xs = x(c0:c0+nw-1); ex = sum(xs.^2); lags = (prev-S):(prev+S); cc = -2*ones(size(lags));
      for q = 1:numel(lags)
        i0 = c0 + lags(q); if i0 < 1 || i0 + nw - 1 > ns, continue; end
        ys = yy(i0:i0+nw-1); cc(q) = sum(xs.*ys)/sqrt(ex*sum(ys.^2));
      end
      [~, q] = max(cc); d = 0;
      if q > 1 && q < numel(lags), den = cc(q-1) - 2*cc(q) + cc(q+1); if den < 0, d = 0.5*(cc(q-1) - cc(q+1))/den; end, end
      prev = lags(q); zcc = (c0 - 1 + nw/2)*dzs; zc(end+1) = zcc; uu(end+1) = (lags(q) + d)*dzs - sh*zcc; %#ok<AGROW>
    end
    ne = numel(zc) - N + 1; ee = zeros(1,ne); ze = zeros(1,ne);
    for i = 1:ne, zz = zc(i:i+N-1); u2 = uu(i:i+N-1); zm = mean(zz); ee(i) = -sum((zz-zm).*(u2-mean(u2)))/sum((zz-zm).^2); ze(i) = zm; end
    sel = ze > 0.004 & ze < 0.036; SN(comp+1, ie) = mean(ee(sel))/std(ee(sel));
  end
  fprintf('strain %6.2f %%: SNRe = %7.1f (plain), %8.1f (companded)\n', e*100, SN(1,ie), SN(2,ie));
end
figure(2); loglog(epsList*100, SN(1,:), 'r-o', epsList*100, SN(2,:), 'g--s'); grid on;
xlabel('applied strain (%)'); ylabel('SNR_e'); legend('plain', 'companded'); title('Strain filter, 2 mm kernel');
