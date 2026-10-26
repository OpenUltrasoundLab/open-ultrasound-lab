%% m33_distortion_matrix_svd.m
% Reflection matrix, focused basis and distortion matrix of a speckle medium
% behind a phase screen; aberration law from the first eigenvector, iteration,
% and the eigenvalue spectrum when the screen is away from the array.
% Purpose : Monochromatic (5 MHz) model, 48 elements, pitch 0.3 mm. Scatterers
%           (complex Gaussian, 400 of them, plus one bright point at x = 0) lie
%           in a thin layer at z = 20 mm, standing in for a time gate. A random
%           delay screen (rms sig_tau, correlation length l) sits at depth d;
%           the straight ray from element x_e to (x_p, z) crosses it at
%           x_e + (x_p - x_e) d/z. Green's function G = exp(-j k r + j phi).
%             K    = G diag(a) G.'                    (element-basis reflection matrix)
%             R_ux = K conj(G0) diag-corrected        (focused on transmit at x_in)
%             D    = R_ux .* conj(G0)                 (remove the geometric wavefront)
%             C    = D D'                             (= U S^2 U', SVD of D)
%             phi_hat = angle(first eigenvector of C), by power iteration.
%           Iteration: correct the input side with phi_hat (post-processing,
%           no new acquisition) and re-estimate. Focusing quality Gamma: the
%           confocal (zero-offset) fraction of the common-midpoint intensity
%           profile of R_xx, normalized by its value without the screen
%           (a Strehl-like ratio; not the FWHM-based F of Lambert et al.).
% Change  : sig_tau, ell, d_scr, Ns, n_iter.
% Expect  : d = 0: residual phase falls from ~1 rad to a few hundredths in
%           3 iterations, Gamma from ~0.5 to 1, the singular spectrum of D
%           collapses onto sigma_1; d = 8 mm (anisoplanatic): a slower singular
%           spectrum, and local windows reach a higher Gamma than one global law.
%           For d > 0 Gamma, not the residual against one point's law, is the
%           metric: a window law may focus well on slightly displaced points.
% Runtime : about 3 s in Octave.
clear; close all; rng(3);

c = 1540; f0 = 5e6; lam = c/f0; k = 2*pi/lam;
Ne = 48; p = 0.3e-3; xe = ((0:Ne-1) - (Ne-1)/2)*p; z = 20e-3;
xin = (-6e-3:0.2e-3:6e-3); Nx = numel(xin);
Ns = 400; sig_tau = 40e-9; ell = 2e-3; n_iter = 4;
xs_sc = [(rand(Ns, 1) - 0.5)*18e-3; 0]; a_sc = [(randn(Ns, 1) + 1j*randn(Ns, 1))/sqrt(2); 6];   % speckle + point at x = 0

dxs = 0.05e-3; xg = -30e-3:dxs:30e-3; gk = exp(-2*(((-3*ell):dxs:(3*ell))/ell).^2); gk = gk/sqrt(sum(gk.^2));
scr = 2*pi*f0*sig_tau*conv(randn(1, numel(xg)), gk, 'same');                    % phase screen phi(x) on the grid (rad)
phiAt = @(x) interp1(xg, scr, x, 'linear');
detr = @(v) v - [ones(numel(v), 1) xe(:)]*([ones(numel(v), 1) xe(:)]\v(:));
resid = @(v) detr(unwrap(angle(exp(1j*v(:)))));                                % wrapped residual -> unwrapped, mean and ramp removed
G0 = exp(-1j*k*sqrt((xe' - xin).^2 + z^2));                                     % Ne x Nx free-space
Ss = @(v) abs(mean(exp(1j*v)))^2;                                               % Strehl of a residual phase vector

for d_scr = [0 8e-3]
    fprintf('\n===== screen at d = %.0f mm (sigma_tau = %.0f ns -> %.2f rad rms, l = %.0f mm) =====\n', d_scr*1e3, sig_tau*1e9, 2*pi*f0*sig_tau, ell*1e3);
    G = exp(-1j*k*sqrt((xe' - xs_sc').^2 + z^2)).*exp(1j*phiAt(xe' + (xs_sc' - xe')*d_scr/z));   % Ne x (Ns+1)
    K = G*diag(a_sc)*G.';                                                        % reflection matrix, element basis
    phiTrue = @(xp) phiAt(xe' + (xp - xe')*d_scr/z);                             % true law for a focal point xp
    % ideal reference width w0 from the screen-free medium
    Gf = exp(-1j*k*sqrt((xe' - xs_sc').^2 + z^2)); Kf = Gf*diag(a_sc)*Gf.';
    dl = -2e-3:0.05e-3:2e-3;                                                     % output offsets for the CMP profile (m)
    cmp = @(Kc, cols) arrayfun(@(dd) mean(abs(sum(conj(exp(-1j*k*sqrt((xe' - xin(cols) - dd).^2 + z^2))).*(Kc*conj(G0(:, cols))), 1)).^2), dl);
    cfr = @(I) I(dl == 0)/sum(I);                                                % confocal fraction of the CMP energy
    w0 = cfr(cmp(Kf, 1:Nx));
    windows = {1:Nx}; Gfin = []; gain = [];                                                            % global window
    if d_scr > 0, windows = [windows, {1:21, 21:41, 41:61}]; end
    for wi = 1:numel(windows)
        cols = windows{wi}; xc = mean(xin(cols));
        psi = zeros(Ne, 1);                                                      % current correction (input and output)
        if wi == 1, fprintf(' window: all %d focal points\n', numel(cols)); else, fprintf(' window %d: x = %.1f to %.1f mm\n', wi - 1, xin(cols(1))*1e3, xin(cols(end))*1e3); end
        for it = 0:n_iter
            Kc = diag(exp(-1j*psi))*K*diag(exp(-1j*psi));                       % correction applied to both sides
            Rux = Kc*conj(G0(:, cols));                                          % focused on transmit
            Dm = Rux.*conj(G0(:, cols));                                         % distortion matrix
            C = Dm*Dm'; Ch = C;                                               % (an entrywise-normalized C./|C| is not positive semidefinite)
            v = ones(Ne, 1)/sqrt(Ne);                                            % power iteration
            for q = 1:60, v = Ch*v; v = v/norm(v); end
            lamC = sort(real(eig(Ch)), 'descend');
            Rxx = G0(:, cols)'*Kc*conj(G0(:, cols)); F = cfr(cmp(Kc, cols))/w0;
            sv = sqrt(max(lamC, 0)); sv = sv/sv(1);                              % singular values of D, normalized
            if d_scr == 0
                err = resid(phiTrue(xc) - psi);
                fprintf('  iter %d: residual phase %.3f rad rms (Strehl %.3f), Gamma = %.2f, sigma2/sigma1 = %.3f\n', it, std(err), Ss(err), F, sv(2));
            elseif it == 0 || it == n_iter
                fprintf('  iter %d: Gamma = %.2f, sigma2/sigma1 = %.3f\n', it, F, sv(2));
            end
            if it == 0 && wi == 1
                pn = lamC/sum(lamC); pn = pn(pn > 0); H = -sum(pn.*log(pn));
                fprintf('   singular values of D (sigma_i/sigma_1): %s ... effective rank exp(H) = %.2f\n', sprintf('%.3f ', sv(1:6)), exp(H));
                conf0 = abs(diag(Rxx)); Rxx0 = Rxx;
            end
            if it < n_iter, psi = psi + detr(unwrap(angle(v))); end               % update (mean and ramp removed: a ramp only shifts the image)
        end
        Gfin(wi) = F; i0 = find(abs(xin(cols)) < 1e-9);
        if ~isempty(i0), cf = abs(diag(Rxx)); gain(wi) = cf(i0); end
        if wi == 1, RxxN = Rxx; g0 = conf0(abs(xin) < 1e-9); end
    end
    fprintf(' bright point at x = 0: confocal amplitude gain, global law %.2f (%.1f dB)\n', gain(1)/g0, 20*log10(gain(1)/g0));
    if d_scr > 0
        fprintf(' bright point at x = 0: confocal amplitude gain, law of window 2 %.2f (%.1f dB)\n', gain(3)/g0, 20*log10(gain(3)/g0));
        fprintf(' final Gamma: global law %.2f, mean of the three local laws %.2f\n', Gfin(1), mean(Gfin(2:4)));
    end
    if d_scr == 0
        figure('Position', [50 50 1100 330]);
        subplot(1, 3, 1); imagesc(xin*1e3, xin*1e3, abs(Rxx0)); axis image; colormap(gray); xlabel('x_{in} (mm)'); ylabel('x_{out} (mm)'); title('|R_{xx}| before');
        subplot(1, 3, 2); imagesc(xin*1e3, xin*1e3, abs(RxxN)); axis image; xlabel('x_{in} (mm)'); title(sprintf('|R_{xx}| after %d iterations', n_iter));
        subplot(1, 3, 3); plot(xe*1e3, detr(phiTrue(0)), 'k', xe*1e3, detr(phiTrue(0)) - resid(phiTrue(0) - psi), 'r--'); xlabel('element (mm)'); ylabel('phase (rad)'); title('aberration law: true vs estimate'); legend('true', 'estimate');
    end
end
