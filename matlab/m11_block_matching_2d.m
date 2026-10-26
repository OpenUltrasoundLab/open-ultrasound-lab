%% m11_block_matching_2d.m
% Purpose : 2D speckle tracking by block matching. An RF speckle image
%           (axial: f0 carrier under a Gaussian envelope; lateral: Gaussian
%           beam, no carrier) is translated by a random subsample shift.
%           A kernel from frame 1 is matched over a search region in
%           frame 2 with 2D NCC and with SSD; the integer peak is refined
%           with (a) separable 1D parabolic fits through the integer peak,
%           (b) axial-first: an axial parabola in each of the 3 lateral
%           columns, then a lateral parabola through the 3 interpolated peak
%           heights, (c) a 2D quadratic least-squares fit (with cross term)
%           on the 3x3 neighbourhood.
% Change  : FWHMx (lateral beam width), dx (line spacing), snrdB,
%           kernel size, useEnv (track the envelope instead of RF).
% Expect  : f0 = 5 MHz, 60 % BW, lateral FWHM 0.6 mm, 20 dB, RF tracking
%           (rms axial / lateral error): separable 0.59 / 15.5 um,
%           axial-first 0.21 / 2.7 um, 2D quadratic LSQ 0.19 / 4.8 um,
%           SSD = NCC. Almost the same at 60 dB SNR: the errors come from
%           the finite speckle kernel and from axial-lateral coupling, not
%           from electronic noise. The separable lateral fit uses the integer
%           axial row, up to half a sample off the RF peak; there the carrier
%           makes the lateral profile of one speckle realization tilted and
%           distorted. Axial-first gives a lateral/axial ratio of about 13,
%           close to the curvature prediction sqrt(R''_ax/R''_lat) = 14.9.
%           Envelope tracking (useEnv = 1) makes the axial error 3-4x worse
%           (no carrier: curvature ratio only 2.6).
rng(1);

%% Parameters
c = 1540;  f0 = 5e6;  fs = 40e6;  B6 = 0.6;    % (m/s), (Hz), (Hz), -6 dB fractional BW
dz = c / (2*fs);                               % axial pixel (m) = 19.25 um
dx = 0.1e-3;                                   % lateral line spacing (m)
FWHMx = 0.6e-3;                                % lateral pulse-echo amplitude FWHM (m)
Nz = 256;  Nx = 64;                            % image size (pixels)
kz = 80;   kx = 15;                            % kernel: 80 axial x 15 lateral pixels (1.54 x 1.5 mm)
sz = 6;    sx = 4;                             % search half-range (pixels)
snrdB = 20;  nTrial = 150;  useEnv = 0;        % SNR, trials, 1 = track |IQ| instead of RF

%% Spectral model (one-sided in axial frequency so that RF and IQ come together)
fz = (0:Nz-1)' / (Nz*dz);  fz(fz >= 1/(2*dz)) = fz(fz >= 1/(2*dz)) - 1/dz;   % axial spatial freq (cycles/m)
fx = (0:Nx-1) / (Nx*dx);   fx(fx >= 1/(2*dx)) = fx(fx >= 1/(2*dx)) - 1/dx;   % lateral spatial freq
fz0 = 2*f0/c;                                  % axial carrier (cycles/m), two-way
sfz = B6*fz0/2/sqrt(2*log(2));                 % axial spectral std (cycles/m)
sigx = FWHMx / (2*sqrt(2*log(2)));             % lateral PSF std (m)
Hz = exp(-(fz - fz0).^2 / (2*sfz^2)) .* (fz > 0);
Hx = exp(-2*pi^2*sigx^2*fx.^2);                % FT of the Gaussian lateral PSF
H  = Hz * Hx;                                  % Nz x Nx
fitPar = @(Rm, R0, Rp) 0.5 * (Rm - Rp) ./ (Rm - 2*R0 + Rp);
% 2D quadratic LSQ on a 3x3 patch: P(u,v) = a + b u + c v + d u^2 + e v^2 + g u v
[uu, vv] = ndgrid(-1:1, -1:1);
Q = [ones(9,1) uu(:) vv(:) uu(:).^2 vv(:).^2 uu(:).*vv(:)];  Qp = (Q'*Q) \ Q';

z0 = 80;  x0 = 25;                             % kernel top-left (1-based)
iz = z0:z0+kz-1;  ix = x0:x0+kx-1;
errSep = zeros(nTrial, 2);  errAF = errSep;  errQ = errSep;  errSSD = errSep;
for t = 1:nTrial
    d = [rand rand];                           % true shift (axial, lateral) in pixels
    W0 = (randn(Nz,Nx) + 1j*randn(Nz,Nx)) .* H;
    a1 = ifft2(W0);
    a2 = ifft2(W0 .* exp(-1j*2*pi*(fz*dz*d(1) + fx*dx*d(2))));
    sc = sqrt(mean(real(a1(:)).^2));
    a1 = a1/sc + (randn(Nz,Nx) + 1j*randn(Nz,Nx)) * 10^(-snrdB/20);   % complex noise: RF part has SNR = snrdB
    a2 = a2/sc + (randn(Nz,Nx) + 1j*randn(Nz,Nx)) * 10^(-snrdB/20);
    if useEnv, I1 = abs(a1); I2 = abs(a2); else, I1 = real(a1); I2 = real(a2); end
    ker = I1(iz, ix);  ker = ker - mean(ker(:));  ek = sum(ker(:).^2);
    R = zeros(2*sz+1, 2*sx+1);  SSD = R;
    for p = -sz:sz
        for q = -sx:sx
            blk = I2(iz+p, ix+q);  b0 = blk - mean(blk(:));
            R(p+sz+1, q+sx+1) = sum(ker(:).*b0(:)) / sqrt(ek * sum(b0(:).^2));
            df = I1(iz, ix) - blk;  SSD(p+sz+1, q+sx+1) = sum(df(:).^2);
        end
    end
    R(:, [1 end]) = -Inf;  R([1 end], :) = -Inf;           % keep the peak off the border
    [~, k] = max(R(:));  [ip, iq] = ind2sub(size(R), k);
    P = R(ip-1:ip+1, iq-1:iq+1);
    sep = [fitPar(P(1,2), P(2,2), P(3,2)), fitPar(P(2,1), P(2,2), P(2,3))];
    % axial-first: axial parabola in each column, lateral parabola through the vertex heights,
    % then the axial offset interpolated to that lateral position (follows a tilted ridge)
    da = fitPar(P(1,:), P(2,:), P(3,:));                     % axial offsets of the 3 columns
    pk = P(2,:) - 0.25 * (P(1,:) - P(3,:)) .* da;            % parabola vertex heights
    u  = fitPar(pk(1), pk(2), pk(3));                        % lateral offset
    af = [da(2) + 0.5*(da(3) - da(1))*u + 0.5*(da(1) - 2*da(2) + da(3))*u^2, u];   % axial offset read at that lateral position
    co = Qp * P(:);                                          % [a b c d e g]
    uv = -[2*co(4) co(6); co(6) 2*co(5)] \ [co(2); co(3)];   % stationary point of the quadratic
    if any(~isfinite(uv)) || any(abs(uv) > 1), uv = sep'; end
    base = [ip - sz - 1, iq - sx - 1];
    errSep(t, :) = base + sep - d;
    errAF(t, :)  = base + af - d;
    errQ(t, :)   = base + uv' - d;
    S2 = -SSD;  S2(:, [1 end]) = -Inf;  S2([1 end], :) = -Inf;
    [~, k] = max(S2(:));  [jp, jq] = ind2sub(size(S2), k);  Ps = -SSD(jp-1:jp+1, jq-1:jq+1);
    errSSD(t, :) = [jp - sz - 1, jq - sx - 1] + [fitPar(Ps(1,2), Ps(2,2), Ps(3,2)), fitPar(Ps(2,1), Ps(2,2), Ps(2,3))] - d;
end

%% Results (um)
pix = [dz dx] * 1e6;
rs = @(e) sqrt(mean(e.^2, 1)) .* pix;        % rms error per direction (um)
modes = {'RF', 'envelope'};
fprintf('pixels: axial %.2f um, lateral %.0f um; lateral FWHM %.2f mm; SNR %d dB; %s tracking\n', pix, FWHMx*1e3, snrdB, modes{useEnv+1});
fprintf('rms error (um)            axial    lateral   ratio\n');
e = rs(errSep); fprintf('NCC + separable parabola  %6.2f   %7.2f   %5.1f\n', e, e(2)/e(1));
e = rs(errAF);  fprintf('NCC + axial-first fits   %6.2f   %7.2f   %5.1f\n', e, e(2)/e(1));
e = rs(errQ);   fprintf('NCC + 2D quadratic LSQ    %6.2f   %7.2f   %5.1f\n', e, e(2)/e(1));
e = rs(errSSD); fprintf('SSD + separable parabola  %6.2f   %7.2f   %5.1f\n', e, e(2)/e(1));
% curvature prediction for RF: -R''(0) axial = (2 pi fz0)^2 + 1/(2 sz^2), lateral = 1/(2 sigx^2)
sigz = 1/(2*pi*sfz);
cAx = (2*pi*fz0)^2 + 1/(2*sigz^2);  cEnv = 1/(2*sigz^2);  cLat = 1/(2*sigx^2);
fprintf('curvature prediction: lateral/axial = %.1f (RF), %.1f (envelope, |IQ| ACF approx.)\n', sqrt(cAx/cLat), sqrt(cEnv/cLat));

%% Figure: last trial
figure('Position', [100 100 1000 420]);
subplot(1, 3, 1);
imagesc(((1:Nx)-1)*dx*1e3, ((1:Nz)-1)*dz*1e3, real(a1)); colormap(gray); hold on;
plot(([x0 x0+kx x0+kx x0 x0]-1)*dx*1e3, ([z0 z0 z0+kz z0+kz z0]-1)*dz*1e3, 'y-', 'LineWidth', 1.5); hold off;
xlabel('lateral (mm)'); ylabel('depth (mm)'); title('RF frame 1 + kernel');
subplot(1, 3, 2);
Rv = R;  Rv(~isfinite(Rv)) = NaN;
imagesc((-sx:sx)*dx*1e3, (-sz:sz)*dz*1e3, Rv); colorbar; hold on;
plot(d(2)*dx*1e3, d(1)*dz*1e3, 'w+', 'MarkerSize', 12, 'LineWidth', 2); hold off;
xlabel('lateral lag (mm)'); ylabel('axial lag (mm)'); title('2D NCC surface');
subplot(1, 3, 3);
plot(errSep(:,2)*pix(2), errSep(:,1)*pix(1), '.', errAF(:,2)*pix(2), errAF(:,1)*pix(1), 'o'); grid on;
xlabel('lateral error (um)'); ylabel('axial error (um)'); title('Error scatter');
legend('separable', 'axial-first', 'Location', 'southoutside');
