%% m35_diffusion.m  --  Perona-Malik and SRAD: flux, stability and when to stop
% Purpose : (1) Flux phi(s) = s g(s) of the two Perona-Malik conductances: maximum at
%           s = kappa (rational g2) and kappa/sqrt(2) (exponential g1), eq. (35.8).
%           (2) Explicit scheme (35.9): with dt <= 1/4 no new extrema appear; dt = 0.3
%           breaks it. (3) On a 128 x 128 speckle phantom (cyst, -6 dB lesion, +10 dB
%           block, point), run Perona-Malik on the log image and SRAD (35.10) on the
%           envelope for 300 iterations and print speckle SNR, edge 10-90% width,
%           lesion gCNR and point -6 dB width at selected iterations.
% Change  : kappa, dt, iteration list, q0 handling.
% Expect  : speckle SNR and gCNR climb fast and saturate; edge width and point width
%           keep growing: the stopping time is a design parameter, not a convergence.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(352);

%% 1. flux maxima
s = linspace(0, 5, 500001); kap = 1;
[~, i1] = max(s.*exp(-(s/kap).^2)); [~, i2] = max(s./(1 + (s/kap).^2));
fprintf('flux maximum: g1 at s = %.4f kappa (1/sqrt2 = %.4f), g2 at s = %.4f kappa\n', s(i1), 1/sqrt(2), s(i2));

%% 2. stability of the explicit scheme on a noisy step (max over the image must not grow)
u0 = [zeros(1, 32) 10*ones(1, 32)] + randn(1, 64);
for dt = [0.25 0.3]
    u = repmat(u0, 64, 1) + randn(64, 64); mx0 = max(u(:)); mn0 = min(u(:));
    for it = 1:200
        Up = u([1 1:63], :) - u; Dn = u([2:64 64], :) - u; Lf = u(:, [1 1:63]) - u; Rt = u(:, [2:64 64]) - u;
        u = u + dt*(Up + Dn + Lf + Rt);                                  % g = 1: the heat equation, worst case
    end
    fprintf('dt = %.2f: range before [%.2f, %.2f], after 200 iterations [%.3g, %.3g]\n', dt, mn0, mx0, min(u(:)), max(u(:)));
end

%% 3. phantom (128 x 128, 0.1 mm, 5 MHz, 60 percent bandwidth, F#2)
N = 128; d = 0.1e-3; c0 = 1540; f0 = 5e6; lam = c0/f0; B = 0.6*f0; sf = B/(2*sqrt(2*log(2))); Fn = 2;
x = ((0:N-1) - N/2 + 0.5)*d; z = ((0:N-1) + 0.5)*d; [X, Z] = meshgrid(x, z);
fk = ([0:N/2-1, -N/2:-1])/(N*d); [KX, KZ] = meshgrid(fk, fk);
H = exp(-(c0*KZ/2).^2/(2*sf^2)) .* max(0, 1 - abs(KX)*lam*Fn);
R = ones(N); R(hypot(X + 3e-3, Z - 3.5e-3) < 1.8e-3) = 0; R(hypot(X - 3e-3, Z - 3.5e-3) < 1.8e-3) = 10^(-0.6); R(X >= 0.5e-3 & Z >= 7.5e-3) = 10;
sc = sqrt(R + 1e-4).*(randn(N) + 1i*randn(N))/sqrt(2); [~, pq] = min(abs(x + 3.5e-3)); [~, pr] = min(abs(z - 9.8e-3));
sc(pr, pq) = sc(pr, pq) + sqrt(1e3*N*N*sum(H(:).^2))/sum(H(:));
A0 = abs(ifft2(fft2(sc).*H))*sqrt(N*N/sum(H(:).^2)); D0 = 20*log10(A0 + 1e-9);
inL = hypot(X - 3e-3, Z - 3.5e-3) < 1.44e-3; outL = hypot(X - 3e-3, Z - 3.5e-3) > 2.34e-3 & hypot(X - 3e-3, Z - 3.5e-3) < 3.42e-3;
bg = X >= -6e-3 & X <= -0.5e-3 & Z >= 6e-3 & Z <= 8.6e-3 & hypot(X + 3.5e-3, Z - 9.8e-3) > 1e-3;
edgesf = @(p) unique(p(round((0:64)*(numel(p) - 1)/64) + 1)); hm = @(h) [h(1:end-2); h(end-1) + h(end)];
gc3 = @(a, b, e) 1 - sum(min(hm(histc(a, e))/numel(a), hm(histc(b, e))/numel(b))); gcnrf = @(a, b) gc3(a(:), b(:), edgesf(sort([a(:); b(:)])));
ec = abs(x - 0.5e-3) <= 2e-3; er = z >= 8.2e-3 & z <= 12.5e-3; xe = (x(ec) - 0.5e-3)'; sg = 0.02e-3*1.03.^(0:130); x0g = (-8:8)*0.03e-3;

%% 4. run both diffusions, metrics every 10 iterations
its = 0:10:300; M = zeros(numel(its), 4, 2); g2 = @(t) 1./(1 + (t/4).^2); dt = 0.2;
for meth = 1:2
    if meth == 1, u = D0; else, u = A0; end
    k = 1;
    for it = 0:300
        if any(it == its)
            if meth == 1, Dk = u; Ak = 10.^(u/20); else, Ak = u; Dk = 20*log10(u + 1e-9); end
            prof = mean(Dk(er, ec), 1)'; best = Inf;
            for is = 1:numel(sg), for ix = 1:numel(x0g)
                Mm = [ones(size(xe)) 0.5*(1 + erf((xe - x0g(ix))/(sg(is)*sqrt(2))))]; e = sum((prof - Mm*(Mm\prof)).^2);
                if e < best, best = e; w = 2.5631*sg(is); end
            end, end
            blk = Ak(pr-3:pr+3, pq-3:pq+3); [~, im] = max(blk(:)); [br, bq] = ind2sub([7 7], im); br = br + pr - 4; bq = bq + pq - 4;
            v = Ak(br, max(1, bq-35):min(N, bq+35)); i0 = bq - max(1, bq-35) + 1; th = v(i0)/2; r = i0; while r < numel(v) && v(r) > th, r = r + 1; end; l = i0; while l > 1 && v(l) > th, l = l - 1; end
            if v(r) > th || v(l) > th, fw = NaN; else, fw = ((r - 1 + (v(r-1) - th)/(v(r-1) - v(r))) - (l + (th - v(l))/(v(l+1) - v(l))))*d; end
            M(k, :, meth) = [mean(Ak(bg))/std(Ak(bg), 1), w*1e3, gcnrf(Ak(inL), Ak(outL)), fw*1e3]; k = k + 1;
        end
        if it == 300, break; end
        if meth == 1
            Up = u([1 1:N-1], :) - u; Dn = u([2:N N], :) - u; Lf = u(:, [1 1:N-1]) - u; Rt = u(:, [2:N N]) - u;
            u = u + dt*(g2(Up).*Up + g2(Dn).*Dn + g2(Lf).*Lf + g2(Rt).*Rt);
        else
            q0s = var(u(bg), 1)/mean(u(bg))^2;
            dN = u([1 1:N-1], :) - u; dS = u([2:N N], :) - u; dW = u(:, [1 1:N-1]) - u; dE = u(:, [2:N N]) - u;
            G2 = (dN.^2 + dS.^2 + dW.^2 + dE.^2)./u.^2; Lp = (dN + dS + dW + dE)./u; q2 = (0.5*G2 - Lp.^2/16)./(1 + 0.25*Lp).^2;
            cq = min(1, max(0, 1./(1 + (q2 - q0s)/(q0s*(1 + q0s))))); u = u + dt*(cq.*dN + cq([2:N N], :).*dS + cq.*dW + cq(:, [2:N N]).*dE);
        end
    end
end
lab = {'Perona-Malik g2, kappa 4 dB/px (log)', 'SRAD, q0 from flat ROI (envelope)'};
for meth = 1:2
    fprintf('\n%s\n  iter  speckle SNR  edge 10-90 (mm)  lesion gCNR  point FWHM (mm)\n', lab{meth});
    for it = [0 20 60 150 300], k = find(its == it); fprintf('  %4d   %8.2f     %8.2f        %7.3f      %8.2f\n', it, M(k, :, meth)); end
end

%% 5. figure: metrics relative to iteration 0
figure(1); clf; tl = {'speckle SNR', 'edge 10-90%', 'lesion gCNR', 'point FWHM'}; sty = {'-', '--', ':', '-.'}; ttl = {'Perona-Malik g2 (log)', 'SRAD (envelope)'};
for meth = 1:2
    subplot(1, 2, meth); hold on;
    for j = 1:4, plot(its, M(:, j, meth)/M(1, j, meth), sty{j}, 'LineWidth', 1.5); end
    hold off; grid on; box on;
    xlabel('iteration (dt = 0.2)'); ylabel('value / value at iteration 0'); title(ttl{meth}); legend(tl, 'Location', 'northwest'); ylim([0 8]);
end
