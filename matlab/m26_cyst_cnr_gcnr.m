%% m26_cyst_cnr_gcnr.m  --  cyst phantom: contrast, CNR (four definitions) and gCNR
% Purpose : simulate a speckle phantom with a circular cyst (5 MHz, 60 % bandwidth, F#2,
%           Hann apodization), then measure image quality as reviewers see it reported.
%   Contrast   C   = 10 log10(<I_in>/<I_out>)  (dB, mean intensities)
%   CNR (ours)     = |mu_in - mu_out| / sqrt(sigma_in^2 + sigma_out^2), on the ENVELOPE
%   variants       : same formula on intensity, on the dB image, and on the dB image
%                    with sigma_out only in the denominator
%   gCNR           = 1 - OVL,  OVL = sum over bins of min(p_in, p_out)
%                    (Rodriguez-Molares et al., IEEE TUFFC 2020). Bins = quantiles of the
%                    pooled data, so the estimate is exactly invariant to strictly monotone
%                    transforms (log, power law, gain).
%   Rayleigh theory for intensity ratio r = <I_out>/<I_in>:  gCNR = r^(-1/(r-1)) (1 - 1/r)
%   CNR_env theory = (sqrt(pi)/2)(1 - q) / sqrt((1 - pi/4)(1 + q^2)), q = 1/sqrt(r)
% Change  : cyst contrasts cset_dB, radius Rc, dynamic ranges DRs, apodization.
% Expect  : CNR_env saturates near 1.91 for an anechoic cyst, CNR on dB images grows with
%           contrast and with clipping; gCNR is identical for log, power and linear maps
%           and can only drop when the display clips (loss of information).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);
c = 1540;  f0 = 5e6;  lam = c/f0;  B = 0.6*f0;  Fn = 2;
Nz = 256;  Nx = 256;  dz = 40e-6;  dx = 50e-6;      % 10.2 mm x 12.8 mm
uz = [0:Nz/2-1, -Nz/2:-1]'/(Nz*dz);  ux = [0:Nx/2-1, -Nx/2:-1]/(Nx*dx);
sf = B/(2*sqrt(2*log(2)));  Hz = exp(-(c*uz/2).^2/(2*sf^2));
t = linspace(-0.5, 0.5, 1001);                       % aperture coordinate (units of width)
coarray = @(apod) interp1(linspace(-1, 1, 2001), conv(apod, apod)/max(conv(apod, apod)), ux*lam*Fn, 'linear', 0);
H_hann = Hz * coarray(0.5 + 0.5*cos(2*pi*t));
H_rect = Hz * coarray(ones(size(t)));
[X, Z] = meshgrid(((0:Nx-1) - Nx/2)*dx, ((0:Nz-1) - Nz/2)*dz);
Rc = 2.5e-3;  rr = sqrt(X.^2 + Z.^2);
inR  = rr < 0.8*Rc;                                  % ROI inside, inset from the edge
outR = rr > 1.3*Rc & rr < 1.9*Rc;                    % background annulus

% ---- metric helpers (anonymous only) ----
cnr   = @(a, b) abs(mean(a) - mean(b)) / sqrt(var(a) + var(b));
qedge = @(p, nb) unique(p(max(1, round(linspace(1, numel(p), nb + 1)))));
fold  = @(h) [h(1:end-2); h(end-1) + h(end)];        % histc puts x == last edge in its own bin
ovl   = @(a, b, e) sum(min(fold(histc(a(:), e))/numel(a), fold(histc(b(:), e))/numel(b)));
gcnrn = @(a, b, nb) 1 - ovl(a, b, qedge(sort([a(:); b(:)]), nb));
gcnr  = @(a, b) gcnrn(a, b, 100);
gth   = @(r) r.^(-1./(r - 1)) .* (1 - 1./r);
cth   = @(r) sqrt(pi)/2*(1 - 1./sqrt(r)) ./ sqrt((1 - pi/4)*(1 + 1./r));
W = (randn(Nz, Nx) + 1i*randn(Nz, Nx))/sqrt(2);     % one scatterer realization for all cases
noise_dB = -50;                                      % thermal noise power re. background speckle
Nn = (randn(Nz, Nx) + 1i*randn(Nz, Nx))/sqrt(2);
Ps = mean(abs(H_hann(:)).^2);                        % speckle power per pixel for unit reflectivity
Nn = Nn * sqrt(Ps * 10^(noise_dB/10));

%% 1. Metrics vs cyst contrast (log display, no clipping)
cset_dB = [-Inf -20 -10 -6 -3];
fprintf('set C | meas C (dB) | CNR env (theory) | CNR int | CNR dB | CNR dB, sigma_out | gCNR (theory at meas C)\n');
for k = 1:numel(cset_dB)
    refl = ones(Nz, Nx);  refl(rr < Rc) = 10^(cset_dB(k)/20);
    A = abs(ifft2(fft2(W .* refl) .* H_hann) + Nn);  I = A.^2;
    L = 10*log10(I / mean(I(outR)));                 % dB image, background at 0 dB
    Cm = 10*log10(mean(I(inR))/mean(I(outR)));  r = 10^(-Cm/10);
    fprintf('%5.0f |   %7.2f   |  %5.3f (%5.3f)  |  %5.3f  | %5.3f  |      %5.3f        |  %5.3f (%5.3f)\n', ...
        cset_dB(k), Cm, cnr(A(inR), A(outR)), cth(r), cnr(I(inR), I(outR)), cnr(L(inR), L(outR)), ...
        abs(mean(L(inR)) - mean(L(outR)))/std(L(outR)), gcnr(A(inR), A(outR)), gth(r));
    if k == 2, A20 = A; end
    if k == 1, A0 = A; end
end
Ar = abs(ifft2(fft2(W .* (rr >= Rc)) .* H_rect) + Nn);
fprintf('anechoic cyst, rect apodization: C = %.1f dB (Hann above: wider mainlobe tails leak more near the edge,\n   rect sidelobes leak more far from it; noise floor %d dB)\n', ...
    10*log10(mean(Ar(inR).^2)/mean(Ar(outR).^2)), noise_dB);

%% 2. Same image, different display maps and dynamic ranges (cyst at -20 dB)
A = A20;  An = A / sqrt(mean(A(outR).^2));           % background rms = 1
DRs = [80 50 40 30 20];
fprintf('\n-20 dB cyst. Display v = T(max(An, 10^(-DR/20))): CNR and gCNR of the DISPLAYED values\n');
fprintf(' DR (dB) | log: CNR  gCNR | power 0.3: CNR  gCNR | linear: CNR  gCNR\n');
for DR = DRs
    Ac = max(An, 10^(-DR/20));
    Dl = 20*log10(Ac);  Dp = Ac.^0.3;  Dn = Ac;
    fprintf('  %4d   |   %5.3f  %5.3f |      %5.3f  %5.3f    |    %5.3f  %5.3f\n', DR, ...
        cnr(Dl(inR), Dl(outR)), gcnr(Dl(inR), Dl(outR)), cnr(Dp(inR), Dp(outR)), gcnr(Dp(inR), Dp(outR)), ...
        cnr(Dn(inR), Dn(outR)), gcnr(Dn(inR), Dn(outR)));
end

%% 3. Estimator sanity check: two halves of the same background must give gCNR ~ 0
half1 = outR & X < 0;  half2 = outR & X >= 0;
mu = abs(ifft2(abs(H_hann).^2));  mu = mu/mu(1,1);  Ac = sum(mu(:).^2)*dx*dz;   % speckle correlation area (m^2)
fprintf('\nsame-statistics regions (%d and %d pixels, but only about %.0f independent speckle cells each):\n', ...
    sum(half1(:)), sum(half2(:)), sum(half1(:))*dx*dz/Ac);
fprintf('   gCNR should be 0; the excess is sampling noise set by the cell count, not by the pixel count\n');
for nb = [10 20 50 100]
    fprintf('   %3d bins: gCNR(same) = %.3f   gCNR(-20 dB cyst) = %.3f\n', nb, gcnrn(A(half1), A(half2), nb), gcnrn(A(inR), A(outR), nb));
end

%% Figures
figure('Position', [60 60 1150 400]);
ax = ((0:Nx-1) - Nx/2)*dx*1e3;  az = ((0:Nz-1) - Nz/2)*dz*1e3;
subplot(1,3,1); imagesc(ax, az, 20*log10(A0/sqrt(mean(A0(outR).^2))), [-50 10]); colormap(gray); axis image; colorbar;
xlabel('lateral (mm)'); ylabel('depth (mm)'); title('anechoic cyst, 60 dB display');
subplot(1,3,2); imagesc(ax, az, 20*log10(max(An, 10^(-30/20))), [-30 10]); axis image; colorbar;
xlabel('lateral (mm)'); title('-20 dB cyst, 40 dB display');
subplot(1,3,3);
Dl = 20*log10(An);  edges = -60:1:15;  ctr = edges(1:end-1) + 0.5;
hi = histc(Dl(inR), edges);  ho = histc(Dl(outR), edges);
plot(ctr, hi(1:end-1)/sum(inR(:)), 'r-', ctr, ho(1:end-1)/sum(outR(:)), 'b-', 'LineWidth', 1.5); hold on;
area(ctr, min(hi(1:end-1)/sum(inR(:)), ho(1:end-1)/sum(outR(:))), 'FaceColor', [0.7 0.7 0.7]);
xlabel('displayed value (dB)'); ylabel('fraction per dB'); grid on;
title(sprintf('-20 dB cyst: inside (red), outside (blue), overlap (gray); gCNR = %.2f', gcnr(A(inR), A(outR))));
