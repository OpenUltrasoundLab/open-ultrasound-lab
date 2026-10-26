%% m31_localization_methods.m
% Purpose : compare the sub-pixel localizers used in ULM (the PALA family)
%           with the Cramer-Rao bound, then show what overlap does to them.
% Model   : one bubble = Gaussian PSF (amplitude std w = FWHM/2.3548,
%           FWHM lambda lateral, 0.75 lambda axial), unknown complex amplitude
%           (random phase), white circular complex noise per pixel.
%           ENR = |a|^2 * pi * wx * wz / (sigma_pix^2 * pix^2) (integrated SNR).
%           Bound (Module 12): sigma_x >= wx / sqrt(ENR), sigma_z >= wz / sqrt(ENR).
% Estimators, all started from the brightest pixel of |IQ|:
%   WA  : amplitude-weighted average of |IQ| in a 5 x 5 box (lambda/4 pixels)
%   GAU : separable 3-point Gaussian fit (parabola through log|IQ|)
%   RS  : radial symmetry centre (Parthasarathy 2012) on a 7 x 7 patch of |IQ|
%   INT : cubic interpolation of |IQ| on a 0.05-pixel grid, then the maximum
%   MLE : maximum of |<y, g(x0, z0)>|^2 over a 9 x 9 window (grid of 0.1 px
%         over +-1.5 px, then 0.01 px; separable templates), the ML estimate for unknown complex amplitude
% Part 2: two equal bubbles at lateral separation s (in FWHM), ENR 30 dB each,
%         random phases: bias and merging of a single-bubble localizer.
% Change  : pixList, enrList, nTrial, sepList.
% Expect  : MLE reaches the bound from 20 dB at both pixel sizes (0.98x and
%           1.04x at 25 dB); at lambda/4 with white pixel noise GAU is 3.4x,
%           WA 2.9x, INT 4.1x and RS 1.8x the bound at 25 dB; at lambda/2
%           GAU and RS come within 1.2-1.3x. WA floors at 11-15 um (pixel
%           locking: its sub-pixel histogram std is 0.19 px instead of 0.29).
%           Overlap: at 0.5-0.75 FWHM 40-63 % of pairs merge into one peak;
%           at 1 FWHM the estimate is pulled 44 um toward the partner;
%           beyond 2 FWHM the bias vanishes.
% Runtime : about 10-20 s in Octave.
rng(1);

%% Parameters (lengths in wavelengths)
lam_um = 308;
fwx = 1.0; fwz = 0.75; wxl = fwx/2.3548; wzl = fwz/2.3548;
pixList = [1/4 1/2];
enrList = [10 15 20 25 30 35 40];
nTrial = 400;
names = {'WA ', 'GAU', 'RS ', 'INT', 'MLE'}; nM = numel(names);
o3 = @(a, b, c) max(-0.5, min(0.5, 0.5*(log(a) - log(c))/min(log(a) - 2*log(b) + log(c), -eps)));
rmse = zeros(numel(pixList), numel(enrList), nM); frac = cell(1, nM);
tic;
for ip = 1:numel(pixList)
    pix = pixList(ip); wx = wxl/pix; wz = wzl/pix;            % PSF std in pixels
    N = 2*ceil(3*wx) + 9; c0 = (N + 1)/2; g = (1:N) - c0;     % patch, centre index c0
    h = 4;                                                    % MLE half-window (pixels)
    for ie = 1:numel(enrList)
        sigp = sqrt(pi*wx*wz/10^(enrList(ie)/10));
        err = zeros(nTrial, nM); fx = zeros(nTrial, nM);
        for k = 1:nTrial
            x0 = rand - 0.5; z0 = rand - 0.5;                 % true offset from the centre pixel
            Y = exp(1i*2*pi*rand)*exp(-(g' - z0).^2/(2*wz^2))*exp(-(g - x0).^2/(2*wx^2)) + sigp/sqrt(2)*(randn(N) + 1i*randn(N));
            A = abs(Y); A([1:c0-4, c0+4:N], :) = 0; A(:, [1:c0-4, c0+4:N]) = 0;   % search the peak within +-3 pixels
            [~, m] = max(A(:)); [i0, j0] = ind2sub([N N], m); A = abs(Y);
            est = zeros(nM, 1);
            % WA (5 x 5 at lambda/4, 3 x 3 at lambda/2)
            q = 2*(pix < 0.3) + 1*(pix >= 0.3); B = A(i0-q:i0+q, j0-q:j0+q); [jj, ii] = meshgrid(-q:q, -q:q);
            est(1) = j0 + sum(B(:).*jj(:))/sum(B(:)) - c0 + 1i*(i0 + sum(B(:).*ii(:))/sum(B(:)) - c0);
            % GAU
            est(2) = j0 + o3(A(i0,j0-1), A(i0,j0), A(i0,j0+1)) - c0 + 1i*(i0 + o3(A(i0-1,j0), A(i0,j0), A(i0+1,j0)) - c0);
            % RS: gradients at the midpoints of a 7 x 7 patch (Roberts cross), weighted LS intersection of the gradient lines
            Pt = A(i0-3:i0+3, j0-3:j0+3);
            du = Pt(1:6, 2:7) - Pt(2:7, 1:6); dv = Pt(1:6, 1:6) - Pt(2:7, 2:7);
            gx = (du - dv)/2; gy = -(du + dv)/2; g2 = gx.^2 + gy.^2;
            [xm, ym] = meshgrid((0:5) - 2.5, (0:5) - 2.5);
            xc0 = sum(g2(:).*xm(:))/sum(g2(:)); yc0 = sum(g2(:).*ym(:))/sum(g2(:));
            wgt = g2./max(0.3, sqrt((xm - xc0).^2 + (ym - yc0).^2));
            gn = sqrt(g2) + eps; nx = -gy./gn; ny = gx./gn; cc = nx.*xm + ny.*ym;
            M2 = [sum(wgt(:).*nx(:).^2) sum(wgt(:).*nx(:).*ny(:)); sum(wgt(:).*nx(:).*ny(:)) sum(wgt(:).*ny(:).^2)];
            v = M2 \ [sum(wgt(:).*nx(:).*cc(:)); sum(wgt(:).*ny(:).*cc(:))]; v = max(-1.5, min(1.5, v));
            est(3) = j0 + v(1) - c0 + 1i*(i0 + v(2) - c0);
            % INT: cubic interpolation on a 0.05-pixel grid within +-1 pixel
            [Xq, Zq] = meshgrid(-1:0.05:1, -1:0.05:1);
            Ai = interp2(-2:2, (-2:2)', A(i0-2:i0+2, j0-2:j0+2), Xq, Zq, 'cubic');
            [~, m] = max(Ai(:)); est(4) = j0 + Xq(m) - c0 + 1i*(i0 + Zq(m) - c0);
            % MLE: |<y,g>|^2 with separable templates, coarse then fine grid
            Yw = Y(i0-h:i0+h, j0-h:j0+h); gg = (-h:h)';
            cx = 0; cz = 0;
            for stage = 1:2
                if stage == 1, cand = (-15:15)*0.1; else, cand = (-10:10)*0.01; end
                Gx = exp(-(gg - (cx + cand)).^2/(2*wx^2)); Gz = exp(-(gg - (cz + cand)).^2/(2*wz^2));
                C = abs(Gz.'*Yw*Gx).^2 ./ (sum(Gz.^2, 1)'*sum(Gx.^2, 1));
                [~, m] = max(C(:)); [a1, a2] = ind2sub(size(C), m); cz = cz + cand(a1); cx = cx + cand(a2);
            end
            est(5) = j0 + cx - c0 + 1i*(i0 + cz - c0);
            err(k,:) = (real(est).' - x0)*pix*lam_um + 1i*(imag(est).' - z0)*pix*lam_um;
            fx(k,:) = real(est).' - round(real(est).');           % sub-pixel position of the estimate
        end
        rmse(ip, ie, :) = sqrt(mean(real(err).^2, 1));
        if ip == 1 && enrList(ie) == 30, for m = 1:nM, frac{m} = fx(:,m); end, end
    end
end
crlb = wxl*lam_um./sqrt(10.^(enrList/10));
fprintf('Lateral RMSE (um) versus ENR, %d trials per point; CRLB sigma_x = w_x/sqrt(ENR)\n', nTrial);
for ip = 1:numel(pixList)
    fprintf('\npixel = lambda/%d (%.0f um)\n   ENR(dB)  CRLB  %s\n', round(1/pixList(ip)), pixList(ip)*lam_um, sprintf('%8s', names{:}));
    for ie = 1:numel(enrList), fprintf('   %4d   %6.1f %s\n', enrList(ie), crlb(ie), sprintf('%8.1f', squeeze(rmse(ip,ie,:)))); end
    fprintf('   ratio to CRLB at 25 dB: %s\n', sprintf('%8.2f', squeeze(rmse(ip, enrList == 25, :))/crlb(enrList == 25)));
end
fprintf('\nPixel locking at lambda/4, 30 dB: std of the sub-pixel part of x-estimates (uniform = %.3f px)\n', sqrt(1/12));
for m = 1:nM, fprintf('   %s  %.3f px\n', names{m}, std(frac{m})); end

%% Part 2: two overlapping bubbles (lambda/4 pixels, ENR 30 dB each)
pix = 1/4; wx = wxl/pix; wz = wzl/pix; N = 41; c0 = 21; g = (1:N) - c0; h = 4;
sepList = [0.5 0.75 1 1.25 1.5 2 3];                       % lateral separation in FWHM
sigp = sqrt(pi*wx*wz/10^(30/10)); nT2 = 300;
fprintf('\nTwo equal bubbles, lateral separation s, ENR 30 dB each; single-bubble MLE started at the brightest pixel\n');
fprintf('   s (FWHM)  s (um)   merged %%   bias toward partner (um)   RMSE (um)\n');
bias2 = zeros(size(sepList)); rm2 = bias2; mg = bias2;
for is = 1:numel(sepList)
    s = sepList(is)*fwx/pix;                                % separation in pixels
    e = zeros(nT2,1); merged = 0;
    for k = 1:nT2
        xa = -s/2 + rand - 0.5; xb = xa + s; za = rand - 0.5; zb = za;
        Y = exp(1i*2*pi*rand)*exp(-(g' - za).^2/(2*wz^2))*exp(-(g - xa).^2/(2*wx^2)) + ...
            exp(1i*2*pi*rand)*exp(-(g' - zb).^2/(2*wz^2))*exp(-(g - xb).^2/(2*wx^2)) + sigp/sqrt(2)*(randn(N) + 1i*randn(N));
        A = abs(Y);
        % count local maxima above half the global maximum (merging test)
        Cm = A(2:end-1, 2:end-1); pk = Cm > 0.5*max(A(:));
        for di = -1:1, for dj = -1:1, if di || dj, pk = pk & Cm >= A((2:end-1)+di, (2:end-1)+dj); end, end, end
        merged = merged + (sum(pk(:)) < 2);
        % localize the peak nearest bubble a
        [ii, jj] = find(pk); ii = ii + 1; jj = jj + 1; [~, q] = min((jj - c0 - xa).^2 + (ii - c0 - za).^2);
        i0 = ii(q); j0 = jj(q); Yw = Y(i0-h:i0+h, j0-h:j0+h); gg = (-h:h)'; cx = 0; cz = 0;
        for stage = 1:2
            if stage == 1, cand = (-15:15)*0.1; else, cand = (-10:10)*0.01; end
            Gx = exp(-(gg - (cx + cand)).^2/(2*wx^2)); Gz = exp(-(gg - (cz + cand)).^2/(2*wz^2));
            C = abs(Gz.'*Yw*Gx).^2 ./ (sum(Gz.^2, 1)'*sum(Gx.^2, 1));
            [~, m] = max(C(:)); [a1, a2] = ind2sub(size(C), m); cz = cz + cand(a1); cx = cx + cand(a2);
        end
        e(k) = (j0 + cx - c0 - xa)*pix*lam_um;              % positive = toward the partner (which is at +x)
    end
    bias2(is) = mean(e); rm2(is) = sqrt(mean(e.^2)); mg(is) = 100*merged/nT2;
    fprintf('   %5.2f    %5.0f    %6.1f        %8.1f                %7.1f\n', sepList(is), sepList(is)*fwx*lam_um, mg(is), bias2(is), rm2(is));
end
fprintf('Runtime %.1f s\n', toc);

%% Figures
figure('Position', [50 50 1100 380]);
subplot(1,3,1); semilogy(enrList, crlb, 'k-', 'LineWidth', 2); hold on;
mk = {'o-', 's-', 'd-', '^-', 'v-'};
for m = 1:nM, semilogy(enrList, squeeze(rmse(1,:,m)), mk{m}); end
xlabel('ENR (dB)'); ylabel('lateral RMSE (\mum)'); title('\lambda/4 pixels'); legend(['CRLB', names], 'Location', 'southwest'); grid on;
subplot(1,3,2); hist([frac{1} frac{2} frac{5}], -0.5:0.05:0.5); legend('WA', 'GAU', 'MLE');
xlabel('sub-pixel part of x estimate (pixels)'); ylabel('count'); title('pixel locking, 30 dB, \lambda/4');
subplot(1,3,3); plot(sepList, bias2, 'o-', sepList, rm2, 's-'); hold on;
plot(sepList, mg, 'k^--'); xlabel('separation (FWHM)'); ylabel('\mum  /  % merged'); legend('bias toward partner', 'RMSE', '% merged'); grid on;
title('two bubbles, 30 dB each');
