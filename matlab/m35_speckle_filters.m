%% m35_speckle_filters.m  --  one speckle phantom, ten classical filters, five numbers
% Purpose : Simulate a 256 x 256 B-mode phantom (k-space filtered complex Gaussian
%           scatterers: 5 MHz, 60% bandwidth, F#2, 0.1 mm pixels) with an anechoic cyst,
%           a -6 dB lesion, a +10 dB block with a vertical edge and point targets. Apply
%           mean, median, Lee, Kuan, Frost, Perona-Malik, SRAD, bilateral, non-local means
%           and self-guided filtering, and report for each: cyst and lesion gCNR, envelope
%           CNR, speckle SNR, log-speckle std, edge 10-90% width (erf fit) and point -6 dB
%           widths. Equations (35.5)-(35.14) of Module 35.
% Change  : window sizes, Cu, K, kappa, iterations, NLM h; f0, bw, Fn of the phantom.
% Expect  : input speckle SNR near 1.91 and log std near 5.57 dB; every filter raises
%           SNR and gCNR, and pays in edge width or point width.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(35);

%% 1. phantom
N = 256; d = 0.1e-3; c0 = 1540; f0 = 5e6; lam = c0/f0; B = 0.6*f0; sf = B/(2*sqrt(2*log(2))); Fn = 2; z0 = 10e-3;
x = ((0:N-1) - N/2 + 0.5)*d;  z = z0 + ((0:N-1) + 0.5)*d;  [X, Z] = meshgrid(x, z);
fk = ([0:N/2-1, -N/2:-1])/(N*d);  [KX, KZ] = meshgrid(fk, fk);
H = exp(-(c0*KZ/2).^2/(2*sf^2)) .* max(0, 1 - abs(KX)*lam*Fn);
R = ones(N);  R(hypot(X + 6e-3, Z - 17e-3) < 2.5e-3) = 0;  R(hypot(X - 5e-3, Z - 17e-3) < 3e-3) = 10^(-0.6);
R(X >= 0 & Z >= 25e-3 & Z <= 33.5e-3) = 10;
s = sqrt(R + 1e-4) .* (randn(N) + 1i*randn(N))/sqrt(2);
apt = sqrt(1e3*N*N*sum(H(:).^2))/sum(H(:));                        % +30 dB point targets
pts = [-8e-3 30e-3; 2.5e-3 13e-3; 3.5e-3 13e-3];
for k = 1:3, [~, iq] = min(abs(x - pts(k,1))); [~, ir] = min(abs(z - pts(k,2))); s(ir, iq) = s(ir, iq) + apt; end
zf = ifft2(fft2(s).*H) * sqrt(N*N/sum(H(:).^2));
A = abs(zf);  Dlog = 20*log10(A + 1e-9);

%% 2. ROIs and metric helpers
inC = hypot(X + 6e-3, Z - 17e-3) < 2.0e-3;  outC = hypot(X + 6e-3, Z - 17e-3) > 3.25e-3 & hypot(X + 6e-3, Z - 17e-3) < 4.75e-3;
inL = hypot(X - 5e-3, Z - 17e-3) < 2.4e-3;  outL = hypot(X - 5e-3, Z - 17e-3) > 3.9e-3 & hypot(X - 5e-3, Z - 17e-3) < 5.7e-3;
bg = abs(X) <= 11e-3 & Z >= 21.5e-3 & Z <= 24.3e-3;
cnrf = @(a, b) abs(mean(a) - mean(b))/sqrt(var(a, 1) + var(b, 1));
edgesf = @(p) unique(p(round((0:64)*(numel(p) - 1)/64) + 1));       % 64 quantile bins of the pooled values (as in the page)
hm = @(h) [h(1:end-2); h(end-1) + h(end)];                            % values equal to the top edge go to the top bin
gc3 = @(a, b, e) 1 - sum(min(hm(histc(a, e))/numel(a), hm(histc(b, e))/numel(b)));
gcnrf = @(a, b) gc3(a(:), b(:), edgesf(sort([a(:); b(:)])));
ecols = abs(x) <= 4e-3;  erows = z >= 26e-3 & z <= 32.5e-3;  xe = x(ecols)';
sg = 0.02e-3*1.03.^(0:120);  x0g = (-10:10)*0.025e-3;               % erf-fit grid: sigma and edge position
rad = 35;  [~, ipq] = min(abs(x - pts(1,1)));  [~, ipr] = min(abs(z - pts(1,2)));

%% 3. filter helpers (replicate padding, box statistics)
padr = @(I, r) I([ones(1, r) 1:N N*ones(1, r)], [ones(1, r) 1:N N*ones(1, r)]);
boxm = @(I, r) conv2(ones(2*r + 1, 1)/(2*r + 1), ones(1, 2*r + 1)/(2*r + 1), padr(I, r), 'valid');
CuE = sqrt(4/pi - 1);  sn = 5.57;
names = {'input', 'mean 5x5 (env)', 'median 5x5 (log)', 'Lee 7x7 (env)', 'Kuan 7x7 (env)', 'Frost 7x7 K=2 (env)', ...
         'Perona-Malik k=4, 40 it (log)', 'SRAD 60 it (env)', 'bilateral 2 px, 6 dB (log)', 'NLM P5 S13 h=1 (log)', 'guided r4 (log)'};
out = cell(1, numel(names));  dom = zeros(1, numel(names));        % dom 0: envelope output, 1: dB output
out{1} = A;
out{2} = boxm(A, 2);
P = padr(Dlog, 2); st = zeros(N, N, 25); k = 0;
for a = -2:2, for b = -2:2, k = k + 1; st(:, :, k) = P(3 + a:N + 2 + a, 3 + b:N + 2 + b); end, end
out{3} = median(st, 3); dom(3) = 1;
m = boxm(A, 3); v = max(0, boxm(A.^2, 3) - m.^2); ci2 = v./m.^2;
Wl = max(0, 1 - CuE^2./ci2); out{4} = m + Wl.*(A - m);
Wk = Wl/(1 + CuE^2);        out{5} = m + Wk.*(A - m);
P = padr(A, 3); num = zeros(N); den = zeros(N);
for a = -3:3, for b = -3:3, w = exp(-2*ci2*hypot(a, b)); num = num + w.*P(4 + a:N + 3 + a, 4 + b:N + 3 + b); den = den + w; end, end
out{6} = num./den;
u = Dlog; g2 = @(t) 1./(1 + (t/4).^2);
for it = 1:40
    Up = u([1 1:N-1], :) - u; Dn = u([2:N N], :) - u; Lf = u(:, [1 1:N-1]) - u; Rt = u(:, [2:N N]) - u;
    u = u + 0.2*(g2(Up).*Up + g2(Dn).*Dn + g2(Lf).*Lf + g2(Rt).*Rt);
end
out{7} = u; dom(7) = 1;
I = A;
for it = 1:60
    q0s = var(I(bg), 1)/mean(I(bg))^2;
    dN = I([1 1:N-1], :) - I; dS = I([2:N N], :) - I; dW = I(:, [1 1:N-1]) - I; dE = I(:, [2:N N]) - I;
    G2 = (dN.^2 + dS.^2 + dW.^2 + dE.^2)./I.^2; Lp = (dN + dS + dW + dE)./I;
    q2 = (0.5*G2 - Lp.^2/16)./(1 + 0.25*Lp).^2;
    cq = min(1, max(0, 1./(1 + (q2 - q0s)/(q0s*(1 + q0s)))));
    cS = cq([2:N N], :); cE = cq(:, [2:N N]);
    I = I + 0.2*(cq.*dN + cS.*dS + cq.*dW + cE.*dE);
end
out{8} = I;
P = padr(Dlog, 4); num = zeros(N); den = zeros(N);
for a = -4:4, for b = -4:4, Q = P(5 + a:N + 4 + a, 5 + b:N + 4 + b); w = exp(-(a^2 + b^2)/8 - (Q - Dlog).^2/72); num = num + w.*Q; den = den + w; end, end
out{9} = num./den; dom(9) = 1;
pr = 2; sr = 6; P = padr(Dlog, sr); acc = zeros(N); wsum = zeros(N); wmax = zeros(N);
dref = mean([mean((Dlog(bg) - Dlog(find(bg) + 8)).^2), mean((Dlog(bg) - Dlog(find(bg) + 8*N)).^2)]);
for a = -sr:sr, for b = -sr:sr
    if a == 0 && b == 0, continue; end
    Q = P(sr + 1 + a:N + sr + a, sr + 1 + b:N + sr + b); dd = boxm((Dlog - Q).^2, pr); w = exp(-dd/dref);
    acc = acc + w.*Q; wsum = wsum + w; wmax = max(wmax, w);
end, end
out{10} = (acc + wmax.*Dlog)./(wsum + wmax); dom(10) = 1;
mu = boxm(Dlog, 4); vv = max(0, boxm(Dlog.^2, 4) - mu.^2); aa = vv./(vv + sn^2); bb = mu - aa.*mu;
out{11} = boxm(aa, 4).*Dlog + boxm(bb, 4); dom(11) = 1;

%% 4. metrics table
fprintf('%-30s %6s %6s %6s %6s %6s %6s %7s %7s %7s\n', 'filter', 'gCNRc', 'gCNRl', 'CNRc', 'CNRl', 'SNR0', 'logsd', 'edge', 'FWHMl', 'FWHMa');
res = zeros(numel(names), 9);
for k = 1:numel(names)
    if dom(k), Dk = out{k}; Ak = 10.^(Dk/20); else, Ak = out{k}; Dk = 20*log10(Ak + 1e-9); end
    prof = mean(Dk(erows, ecols), 1)';  best = Inf;
    for is = 1:numel(sg), for ix = 1:numel(x0g)
        M = [ones(size(xe)) 0.5*(1 + erf((xe - x0g(ix))/(sg(is)*sqrt(2))))];  e = sum((prof - M*(M\prof)).^2);
        if e < best, best = e; w1090 = 2.5631*sg(is); end
    end, end
    blk = Ak(ipr-3:ipr+3, ipq-3:ipq+3); [~, im] = max(blk(:)); [br, bq] = ind2sub([7 7], im); br = br + ipr - 4; bq = bq + ipq - 4;
    wd = zeros(1, 2);
    for dirn = 1:2
        if dirn == 1, v = Ak(br, bq-rad:bq+rad); else, v = Ak(br-rad:br+rad, bq)'; end
        i0 = rad + 1; th = v(i0)/2; r = i0; while r < numel(v) && v(r) > th, r = r + 1; end; l = i0; while l > 1 && v(l) > th, l = l - 1; end
        if v(r) > th || v(l) > th, wd(dirn) = NaN; else, wd(dirn) = ((r - 1 + (v(r-1) - th)/(v(r-1) - v(r))) - (l + (th - v(l))/(v(l+1) - v(l))))*d; end
    end
    res(k, :) = [gcnrf(Ak(inC), Ak(outC)), gcnrf(Ak(inL), Ak(outL)), cnrf(Ak(inC), Ak(outC)), cnrf(Ak(inL), Ak(outL)), ...
                 mean(Ak(bg))/std(Ak(bg), 1), std(Dk(bg), 1), w1090*1e3, wd*1e3];
    fprintf('%-30s %6.3f %6.3f %6.2f %6.2f %6.2f %6.2f %7.2f %7.2f %7.2f\n', names{k}, res(k, :));
end
fprintf('theory: speckle SNR 1.913, log-speckle std 5.57 dB; widths in mm (edge 10-90%%, point -6 dB amplitude)\n');

%% 5. figure
figure(1); clf; sel = [1 4 7 8 10 11];
for k = 1:6
    j = sel(k); if dom(j), Dk = out{j}; else, Dk = 20*log10(out{j} + 1e-9); end
    subplot(2, 3, k); imagesc(x*1e3, z*1e3, Dk, [-45 15]); axis image; colormap(gray);
    title(names{j}); xlabel('x (mm)'); ylabel('z (mm)');
end
