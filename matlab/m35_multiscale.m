%% m35_multiscale.m  --  speckle cells vs pixels, log speckle, pyramid halos, CLAHE slope bound
% Purpose : (1) Window averaging of correlated speckle, eq. (35.5): Monte Carlo speckle-SNR
%           gain of a w x w mean filter on intensity vs the Gaussian-PSF formula and the
%           cell-count asymptote sqrt(|W|/A_c). (2) Log speckle of L-look averages, eq. (35.4):
%           std 4.343 sqrt(psi1(L)) and mean bias. (3) Unsharp masking of a step, eq. (35.17):
%           overshoot lambda/2. (4) Laplacian-pyramid gains (35.16): overshoot of a blurred step
%           for several gain sets. (5) CLAHE clip (35.18): maximum mapping slope vs beta.
% Change  : speckle cell size (sx, sz), window list, L list, pyramid gains, beta list.
% Expect  : 5 x 5 window gain 1.85 (not 5); log std 5.57/3.49/2.31/1.58/1.10 dB for
%           L = 1/2/4/8/16; unsharp overshoot 0.5 lambda; pyramid edge boost (1, 1.8, 2, 1.4)
%           overshoots by about 22%; CLAHE slope <= (L-1)/B (beta + E/N), the bound reached
%           whenever the peak bin is clipped; no clipping above the peak bin (plain equalization).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(353);

%% 1. speckle cells vs pixels (intensity speckle, Gaussian PSF, 0.1 mm grid)
sx = 0.514/2.3548/0.1; sz = 0.226/2.3548/0.1;            % std of |mu|^2 in pixels (lateral, axial), Module 26
N = 512; [KX, KZ] = meshgrid(([0:N/2-1, -N/2:-1])/N);
H = exp(-(2*pi^2)*(sx^2*KX.^2 + sz^2*KZ.^2));            % |H| whose |mu|^2 is Gaussian with std sx, sz
I = abs(ifft2(fft2(randn(N) + 1i*randn(N)).*H)).^2; I = I/mean(I(:));
f = @(a, s) 2/a^2*(a*s*sqrt(pi/2)*erf(a/(s*sqrt(2))) - s^2*(1 - exp(-a^2/(2*s^2))));
fprintf('single-pixel intensity speckle SNR %.3f (theory 1)\n', mean(I(:))/std(I(:), 1));
fprintf(' w   MC gain  formula  sqrt(|W|/Ac)  cells\n');
wMC = [3 5 7 11 21]; gMC = zeros(size(wMC));
for w = wMC
    Ib = conv2(I, ones(w)/w^2, 'valid'); g = (mean(Ib(:))/std(Ib(:), 1))/(mean(I(:))/std(I(:), 1)); gMC(wMC == w) = g;
    fprintf('%2d   %6.2f   %6.2f   %8.2f     %5.2f\n', w, g, 1/sqrt(f(w, sx)*f(w, sz)), w/sqrt(2*pi*sx*sz), w^2/(2*pi*sx*sz));
end

%% 2. log speckle of L-look averages
fprintf('\n L   log std MC  theory   mean bias MC  theory (dB)\n');
tri = @(L) sum(1./((L:L+1e6-1).^2));                       % trigamma by direct summation
for L = [1 2 4 8 16]
    S = mean(-log(rand(2e5, L)), 2); D = 10*log10(S);
    fprintf('%2d   %8.3f   %6.3f   %10.3f   %7.3f\n', L, std(D), 4.3429*sqrt(tri(L)), mean(D), 4.3429*(psi(L) - log(L)));
end

%% 3. unsharp masking of a step: overshoot lambda/2
xx = (-400:400)*0.01; sig = 1; Gk = exp(-xx.^2/(2*sig^2)); Gk = Gk/sum(Gk); u = double(xx > 0);
ub = conv([zeros(1, 400) u ones(1, 400)], Gk, 'same'); ub = ub(401:1201);
for lamU = [0.5 1 2]
    e = u + lamU*(u - ub); fprintf('unsharp lambda %.1f: overshoot %.3f (lambda/2 = %.3f), value at x = sigma: %.3f (0.159 lambda = %.3f)\n', ...
        lamU, max(e) - 1, lamU/2, e(xx == 1) - 1, 0.1587*lamU);
end

%% 4. Laplacian pyramid (a = 0.4) on a blurred step
w5 = [0.05 0.25 0.4 0.25 0.05]; M = 256; xs = (0:M-1) - M/2 + 0.5; st = 10*0.5*(1 + erf(xs/(2.3*sqrt(2))));
red = @(g) conv(g([1 1 1:end end end]), w5, 'valid'); red2 = @(g) g(1:2:end);
ctr = abs(xs) < 40;                                        % central region: the replicate-padded expand is not exact near the ends
fprintf('\n gains g0..g3        overshoot  undershoot (%% of step, central 8 mm)\n');
for G = [1 1 1 1; 0.2 0.5 1 1; 1 1.8 2 1.4; 1 1 2 1; 3 1 1 1]'
    Gs = cell(1, 5); Ls = cell(1, 4); Gs{1} = st;
    for k = 1:4
        Gs{k+1} = red2(red(Gs{k})); up = zeros(1, 2*numel(Gs{k+1})); up(1:2:end) = Gs{k+1};
        Ls{k} = Gs{k} - 2*conv(up([1 1 1:end end end]), w5, 'valid');
    end
    rec = Gs{5};
    for k = 4:-1:1, up = zeros(1, 2*numel(rec)); up(1:2:end) = rec; rec = 2*conv(up([1 1 1:end end end]), w5, 'valid') + G(k)*Ls{k}; end
    fprintf(' %.1f %.1f %.1f %.1f          %6.1f     %6.1f\n', G, 100*max(0, max(rec(ctr)) - 10)/10, 100*min(0, min(rec(ctr)))/10);
    if all(G == 1), fprintf('   reconstruction error at unit gains, central region: %.1e\n', max(abs(rec(ctr) - st(ctr)))); end
    if all(G' == [1 1.8 2 1.4]), recE = rec; end
end

%% 5. CLAHE clip limit bounds the mapping slope
Dl = 10*log10(-log(rand(4096, 1)));  g = min(255, max(0, round(255*(Dl + 30)/40)));   % one tile of log speckle, 256 levels
fprintf('\n beta   max slope   bound (L-1)/B (beta + E/N)   fraction clipped E/N\n');
for beta = [1 2 4 8 16]
    h = histc(g, 0:255); lim = beta*numel(g)/256; E = sum(max(0, h - lim)); h = min(h, lim) + E/256;
    fprintf(' %5.3g    %6.2f        %6.2f                    %5.3f\n', beta, max(255*h/numel(g)), 255/256*(beta + E/numel(g)), E/numel(g));
end

%% 6. figure
figure(1); clf;
subplot(1, 2, 1); ws = 1:25; plot(ws, arrayfun(@(a) 1/sqrt(f(a, sx)*f(a, sz)), ws), ws, ws/sqrt(2*pi*sx*sz), '--', wMC, gMC, 'ko', 'LineWidth', 1.5); grid on;
xlabel('mean-filter window side w (pixels, 0.1 mm)'); ylabel('speckle SNR gain'); legend('Gaussian-PSF formula', 'sqrt(|W|/A_c)', 'Monte Carlo', 'Location', 'northwest'); title('count speckle cells, not pixels');
subplot(1, 2, 2); plot(xs*0.1, st, '--', xs*0.1, recE, 'LineWidth', 1.5); grid on; xlim([-3 3]);
xlabel('x (mm)'); ylabel('dB'); legend('step', 'pyramid gains 1, 1.8, 2, 1.4', 'Location', 'southeast'); title('pyramid edge boost: halo');
