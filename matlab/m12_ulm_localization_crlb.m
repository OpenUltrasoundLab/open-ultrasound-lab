%% m12_ulm_localization_crlb.m
% Purpose : ULM localization precision for one isolated bubble in a 2D
%           beamformed IQ image, compared with the Cramer-Rao bound.
%   Model : s(x,z) = a * gx(x - x0) * gz(z - z0), Gaussian PSF with
%           lateral/axial FWHM, a complex and unknown (random phase).
%           Noise: circular complex Gaussian, white inside the system band
%           |kx| <= 0.483/wx, |kz| <= 0.483/wz (PSF energy spectrum -40 dB),
%           zero outside, density N0. This mimics beamformed noise, which
%           is band-limited: pixels finer than Nyquist add no information.
%   Bound : sigma_x >= wx / sqrt(ENR), sigma_z >= wz / sqrt(ENR),
%           w = FWHM / 2.3548, ENR = ||s||^2 / N0 = |a|^2 pi wx wz / N0
%           (integrated SNR; complex noise; amplitude and phase unknown).
%   Estimators (all from the pixel with the largest |IQ|):
%           WA  : amplitude-weighted centroid of |IQ| in a +-FWHM box
%           IC  : intensity-weighted centroid of |IQ|^2 in the same box
%           PAR : separable 3-point parabola on |IQ|
%           GAU : separable 3-point Gaussian fit (parabola on log|IQ|)
%           MLE : maximize |<y, g(x,z)>|^2 / ||g||^2 (Newton, 2D)
% Change  : fwhmX, fwhmZ, pixList, enrList, nTrial, f0.
% Expect  : MLE within 0.95-1.05x of the bound for ENR >= 20 dB at every
%           pixel size; IC (intensity centroid) 1.0-1.3x; GAU 1.2-1.4x at
%           lambda/2 but ~2.2x at lambda/4 and lambda/8 (a 3-point fit uses
%           less of the PSF as pixels shrink); PAR has a bias floor at
%           lambda/2 (5.7x at 40 dB); WA 1.2-3x; NS (pixel centre) floors at
%           pix/sqrt(12). Peak-pixel SNR = ENR - 4.8 dB (lambda/4, lambda/8).
%           At 16 dB 1% outliers (a noise peak beats the bubble) already give
%           2-2.5x; at 10 dB about half the trials are outliers (out%).
rng(1);
% 3-point vertex offset in pixels, clamped to +-0.5 (exact bound when the
% centre sample is the maximum; guards against edge peaks at low SNR)
off3 = @(v) max(-0.5, min(0.5, 0.5*(v(1) - v(3))/min(v(1) - 2*v(2) + v(3), -eps)));

%% Parameters (lengths in wavelengths)
f0 = 5e6; c = 1540; lam = c/f0;            % 5 MHz: lambda = 308 um
fwhmX = 1.0; fwhmZ = 0.75;                 % PSF FWHM (lambda), lateral and axial
wx = fwhmX/(2*sqrt(2*log(2))); wz = fwhmZ/(2*sqrt(2*log(2)));
pixList = [1/2 1/4 1/8];                   % square pixel size (lambda)
enrList = [10 16 20 25 30 40];            % integrated SNR (dB)
nTrial = 300;
names = {'MLE', 'IC ', 'WA ', 'GAU', 'PAR', 'NS '};
nEst = numel(names);
crlbX = wx ./ sqrt(10.^(enrList/10));
crlbZ = wz ./ sqrt(10.^(enrList/10));
rmseX = zeros(numel(pixList), numel(enrList), nEst);
rmseZmle = zeros(numel(pixList), numel(enrList));
outMLE = zeros(numel(pixList), numel(enrList));   % fraction of MLE outliers
snrPkOff = zeros(1, numel(pixList));             % ENR - peak-pixel SNR (dB)
tic;
for ip = 1:numel(pixList)
    pix = pixList(ip);
    nx = ceil(3*fwhmX/pix); nz = ceil(3*fwhmZ/pix);
    xs = (-nx:nx)*pix; zs = (-nz:nz)'*pix;            % pixel grid
    Nx = numel(xs); Nz = numel(zs);
    Kx = ceil(fwhmX/pix); Kz = ceil(fwhmZ/pix);       % centroid half-box (pixels)
    Lx = 2^nextpow2(2*Nx); Lz = 2^nextpow2(2*Nz);
    kx = [0:Lx/2, -Lx/2+1:-1]/(Lx*pix); kz = [0:Lz/2, -Lz/2+1:-1]'/(Lz*pix);
    mask = (abs(kz) <= 0.483/wz) * (abs(kx) <= 0.483/wx);   % rectangular band
    Eg = pi*wx*wz;                                    % integral of |g|^2 (a = 1)
    se = zeros(numel(enrList), nEst); seZ = zeros(numel(enrList), 1); nOut = zeros(numel(enrList), 1);
    sig2pix = mean(mask(:))/pix^2;                    % per-pixel noise variance (N0 = 1)
    snrPkOff(ip) = 10*log10(Eg*sig2pix);              % ENR = peak SNR * Eg * sig2pix
    for t = 1:nTrial
        x0 = (rand - 0.5)*pix; z0 = (rand - 0.5)*pix; ph = 2*pi*rand;
        Wn = (randn(Lz, Lx) + 1j*randn(Lz, Lx)) * sqrt(1/(2*pix*pix));   % density N0 = 1
        Wn = ifft2(fft2(Wn) .* mask);
        Wn = Wn(1:Nz, 1:Nx);
        G = exp(-(zs - z0).^2/(2*wz^2)) * exp(-(xs - x0).^2/(2*wx^2));
        for ie = 1:numel(enrList)
            A = sqrt(10^(enrList(ie)/10)/Eg);
            Y = A*exp(1j*ph)*G + Wn;
            M = abs(Y);
            [~, imax] = max(M(:)); [iz, ix] = ind2sub(size(M), imax);
            iz = min(max(iz, 2), Nz - 1); ix = min(max(ix, 2), Nx - 1);
            rz = max(iz - Kz, 1):min(iz + Kz, Nz); rx = max(ix - Kx, 1):min(ix + Kx, Nx);
            B = M(rz, rx); XX = ones(numel(rz), 1)*xs(rx); ZZ = zs(rz)*ones(1, numel(rx));
            wa = [sum(B(:).*XX(:)), sum(B(:).*ZZ(:))]/sum(B(:));
            ic = [sum(B(:).^2.*XX(:)), sum(B(:).^2.*ZZ(:))]/sum(B(:).^2);
            mx = M(iz, ix-1:ix+1); mz = M(iz-1:iz+1, ix).';
            par = [xs(ix) + pix*off3(mx), zs(iz) + pix*off3(mz)];
            gau = [xs(ix) + pix*off3(log(mx)), zs(iz) + pix*off3(log(mz))];
            % MLE: Newton on L(x,z) = ln|gz' Y gx|^2 - ln(||gx||^2 ||gz||^2), start at IC
            p = ic; p0 = p;
            for it = 1:8
                dx = xs.' - p(1); dz = zs - p(2);
                gx = exp(-dx.^2/(2*wx^2)); gx1 = gx.*dx/wx^2; gx2 = gx.*(dx.^2/wx^4 - 1/wx^2);
                gz = exp(-dz.^2/(2*wz^2)); gz1 = gz.*dz/wz^2; gz2 = gz.*(dz.^2/wz^4 - 1/wz^2);
                YX = Y*[gx gx1 gx2];                       % Nz x 3
                Z = gz.'*YX(:,1); Zx = gz.'*YX(:,2); Zxx = gz.'*YX(:,3);
                Zz = gz1.'*YX(:,1); Zzz = gz2.'*YX(:,1); Zxz = gz1.'*YX(:,2);
                P = abs(Z)^2;
                Px = 2*real(conj(Z)*Zx); Pz = 2*real(conj(Z)*Zz);
                Pxx = 2*(abs(Zx)^2 + real(conj(Z)*Zxx)); Pzz = 2*(abs(Zz)^2 + real(conj(Z)*Zzz));
                Pxz = 2*real(conj(Zx)*Zz + conj(Z)*Zxz);
                ex = sum(gx.^2); ex1 = 2*sum(gx.*gx1); ex2 = 2*sum(gx1.^2 + gx.*gx2);
                ez = sum(gz.^2); ez1 = 2*sum(gz.*gz1); ez2 = 2*sum(gz1.^2 + gz.*gz2);
                grad = [Px/P - ex1/ex; Pz/P - ez1/ez];
                H = [Pxx/P - (Px/P)^2 - ex2/ex + (ex1/ex)^2, Pxz/P - Px*Pz/P^2; ...
                     Pxz/P - Px*Pz/P^2, Pzz/P - (Pz/P)^2 - ez2/ez + (ez1/ez)^2];
                if H(1,1) >= 0 || det(H) <= 0, break; end
                st = -(H\grad).';
                if norm(st) > pix, st = st*pix/norm(st); end
                p = p + st;
                if norm(p - p0) > 1.5*pix, p = p0 + (p - p0)*1.5*pix/norm(p - p0); break; end
                if norm(st) < 1e-9, break; end
            end
            est = [p; ic; wa; gau; par; xs(ix), zs(iz)];
            se(ie, :) = se(ie, :) + (est(:, 1).' - x0).^2;
            seZ(ie) = seZ(ie) + (p(2) - z0)^2;
            nOut(ie) = nOut(ie) + (abs(p(1) - x0) > fwhmX/2 || abs(p(2) - z0) > fwhmZ/2);
            if ip == 3 && t == nTrial && enrList(ie) == 20   % keep one example for the figure
                Yex = Y; xsEx = xs; zsEx = zs; x0Ex = x0; z0Ex = z0; estEx = est;
            end
        end
    end
    rmseX(ip, :, :) = reshape(sqrt(se/nTrial), [1, numel(enrList), nEst]);
    rmseZmle(ip, :) = sqrt(seZ.'/nTrial);
    outMLE(ip, :) = nOut.'/nTrial;
end
fprintf('PSF FWHM lateral %.2f lambda, axial %.2f lambda; f0 = %.1f MHz (lambda = %.0f um); %d trials\n', ...
    fwhmX, fwhmZ, f0/1e6, lam*1e6, nTrial);
fprintf('Lateral CRLB sigma_x = wx/sqrt(ENR); entries are RMSE / CRLB (lateral)\n');
for ip = 1:numel(pixList)
    fprintf('\npixel = lambda/%d (%.2f FWHM_x); peak-pixel SNR = ENR - %.1f dB\n', round(1/pixList(ip)), pixList(ip)/fwhmX, snrPkOff(ip));
    fprintf('  ENR(dB)  CRLB_x(lambda)  CRLB_x(um) ');
    fprintf('%7s', names{:}); fprintf('   MLE_z  out%%\n');
    for ie = 1:numel(enrList)
        fprintf('  %6d  %10.4f  %10.2f  ', enrList(ie), crlbX(ie), crlbX(ie)*lam*1e6);
        fprintf('%7.2f', squeeze(rmseX(ip, ie, :)).'/crlbX(ie));
        fprintf('  %6.2f  %4.1f\n', rmseZmle(ip, ie)/crlbZ(ie), 100*outMLE(ip, ie));
    end
end
fprintf('\nFWHM-lambda rule: ENR 20 dB -> sigma_x = lambda/%.1f; 30 dB -> lambda/%.1f\n', ...
    1/(wx/sqrt(100)), 1/(wx/sqrt(1000)));
fprintf('Elapsed %.1f s\n', toc);

%% Figures
figure('Position', [60 60 1150 360]);
subplot(1, 3, 1);
imagesc(xsEx, zsEx, 20*log10(abs(Yex)/max(abs(Yex(:))))); caxis([-30 0]); colormap(gray); axis image; colorbar;
hold on; plot(x0Ex, z0Ex, 'g+', 'MarkerSize', 12, 'LineWidth', 1.5);
plot(estEx(1,1), estEx(1,2), 'ro', estEx(4,1), estEx(4,2), 'cs', 'MarkerSize', 8);
xlabel('lateral x (\lambda)'); ylabel('axial z (\lambda)');
title('|IQ| dB, \lambda/8, ENR 20 dB: + true, o MLE, [] Gauss');
subplot(1, 3, 2);
ip = 2; sty = {'b-o', 'g-s', '-^', 'm-d', 'r-v', 'k:'};
for k = 1:nEst
    semilogy(enrList, squeeze(rmseX(ip, :, k)), sty{k}); hold on;
end
semilogy(enrList, crlbX, 'k-', 'LineWidth', 2);
xlabel('ENR (dB)'); ylabel('lateral RMSE (\lambda)'); grid on;
title('\lambda/4 pixels'); legend([names, {'CRLB'}], 'Location', 'northeast');
subplot(1, 3, 3);
ie = find(enrList == 30);
bar(squeeze(rmseX(:, ie, 1:5))/crlbX(ie));
set(gca, 'XTickLabel', {'\lambda/2', '\lambda/4', '\lambda/8'}); ylim([0 4]);
hold on; plot([0.5 3.5], [1 1], 'k--');
ylabel('RMSE / CRLB (lateral)'); title('ENR 30 dB'); legend(names(1:5), 'Location', 'northeast');
