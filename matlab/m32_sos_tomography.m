%% m32_sos_tomography.m
% Limited-angle echo-mode SoS tomography: two acquisition designs, the
% resolution matrix versus depth, and the singular spectrum.
% Purpose : Unknown: slowness deviation ds = 1/c - 1/c0 on 16 x 14 cells of
%           2.4 mm (x = +-19.2 mm, z = 2 to 35.6 mm) under a 38.4 mm array.
%           Data at every cell centre, built from straight-ray delay errors
%           L_phi(r) (line integral of ds from r to the array at angle phi):
%             'tx-only'  transmit steered, receive fixed at 0 deg; adjacent
%                        transmit angles differenced: d = L_phi(k) - L_phi(k-1)
%             'CMA'      transmit phi_i and receive phi_j from the same angle set;
%                        pairs grouped by mid-angle (phi_i + phi_j)/2 and
%                        adjacent pairs in a group differenced
%           Tikhonov with a gradient penalty; oracle lambda (minimum RMSE
%           against the truth), the best case for each design. R = (A'A + lam^2 G'G)^-1 A'A: row k is the point-spread
%           function of the reconstruction at cell k.
%           Data come from 0.05 mm ray steps through a smooth truth; the
%           inversion uses 0.3 mm steps with nearest-cell weights.
% Change  : thmax, K, noise_ns, the truth.
% Expect  : CMA gives about 2.2x more rows, yet in this straight-ray model
%           every CMA row (W_i + W_j) - (W_i-1 + W_j+1) is a combination of
%           the tx-only rows of the same pixel (printed residual ~1e-15): CMA
%           adds no ray information. Its real-data advantage (higher speckle
%           correlation, no echo-position bias) is outside this model.
%           The error budget is printed per design: with data consistent with
%           the cell model the error comes from noise + regularization; with
%           fine-ray data it is dominated by model (discretization) error,
%           which is larger than the 2 ns noise effect. Both designs smear
%           axially more than laterally and worse with depth; the 5 mm fat
%           layer comes back near 1506 m/s instead of 1450.
% Runtime : about 6 s in Octave.
clear; close all; rng(2);

c0 = 1540; Nx = 16; Nz = 14; h = 2.4;        % reference speed, grid, cell size (mm)
xc = ((0:Nx-1) - (Nx-1)/2)*h; zc = 2 + ((0:Nz-1) + 0.5)*h;
np = Nx*Nz; Xh = 19.2;                       % unknowns, array half-width (mm)
sg = @(r, r0) 1./(1 + exp((r - r0)/0.4));
truth = @(x, z) 1540 + 40*sg(hypot(x - 2, z - 20), 4.5) - 90*(z < 5);   % inclusion + 5 mm fat (1450)
dsf = @(c) (1./c - 1/c0)*1e6;                % ns/mm
[XX, ZZ] = meshgrid(xc, zc); ct = reshape(truth(XX, ZZ)', [], 1);  % index = iz*Nx + ix + 1
cellof = @(x, z) min(max(floor((z - 2)/h), 0), Nz-1)*Nx + min(max(floor((x + Xh)/h), 0), Nx-1) + 1;

% straight-ray integral from point (x0,z0) to the surface at angle phi
% (weights over cells for the inversion; fine numerical value for the data)
ray_w = @(x0, z0, phi, st) accumarray(cellof(x0 - (0:st:z0/cos(phi))'*sin(phi), z0 - (0:st:z0/cos(phi))'*cos(phi)), st, [np 1])';
ray_d = @(x0, z0, phi) sum(dsf(truth(x0 - (0:0.05:z0/cos(phi))*sin(phi), z0 - (0:0.05:z0/cos(phi))*cos(phi))))*0.05;

% gradient penalty
I = speye(np); G = [];
for iz = 0:Nz-1, for ix = 0:Nx-2, k = iz*Nx + ix + 1; G = [G; I(k+1, :) - I(k, :)]; end, end %#ok<AGROW>
for iz = 0:Nz-2, for ix = 0:Nx-1, k = iz*Nx + ix + 1; G = [G; I(k+Nx, :) - I(k, :)]; end, end %#ok<AGROW>
GG = full(G'*G);

thmax = 25; K = 7; noise_ns = 2; phis = linspace(-thmax, thmax, K)*pi/180;
designs = {'tx-only', 'CMA'}; fw = @(p, dd) dd*0.01*sum(interp1(1:numel(p), p, 1:0.01:numel(p)) >= max(p)/2);
fprintf('theta_max = %d deg, K = %d angles, delay noise %.1f ns\n', thmax, K, noise_ns);
fprintf('  design    rows   lambda    trace(R)  RMSE (m/s)  incl. mean (true %.0f)  fat mean (true 1450)   axial / lateral FWHM (mm) at z = 9, 19, 31 mm\n', truth(2, 20));
res = struct(); rowres = 0;
for dsg = 1:2
    Arows = []; d = [];
    for iz = 1:Nz
        for ix = 1:Nx
            x0 = xc(ix); z0 = zc(iz);
            ok = abs(x0 - z0*tan(phis)) <= Xh;            % ray must reach the aperture
            Wp = zeros(K, np); Dp = zeros(K, 1);
            for a = 1:K
                if ok(a), Wp(a, :) = ray_w(x0, z0, phis(a), 0.3); Dp(a) = ray_d(x0, z0, phis(a)); end
            end
            if dsg == 1                                     % tx-only: receive fixed at 0 deg (cancels)
                for a = 2:K
                    if ok(a) && ok(a-1), Arows = [Arows; Wp(a, :) - Wp(a-1, :)]; d = [d; Dp(a) - Dp(a-1)]; end %#ok<AGROW>
                end
            else                                            % CMA: tx i, rx j, group by mid-angle
                Tp = []; for a = 2:K, if ok(a) && ok(a-1), Tp = [Tp; Wp(a, :) - Wp(a-1, :)]; end, end %#ok<AGROW>
                Cp = [];
                for gsum = 2:2*K
                    pr = [];
                    for i = max(1, gsum-K):min(K, gsum-1), j = gsum - i; if i <= j && ok(i) && ok(j), pr = [pr; i j]; end, end %#ok<AGROW>
                    for q = 2:size(pr, 1)
                        w1 = Wp(pr(q, 1), :) + Wp(pr(q, 2), :); w0 = Wp(pr(q-1, 1), :) + Wp(pr(q-1, 2), :);
                        Arows = [Arows; w1 - w0]; d = [d; Dp(pr(q, 1)) + Dp(pr(q, 2)) - Dp(pr(q-1, 1)) - Dp(pr(q-1, 2))]; %#ok<AGROW>
                        Cp = [Cp; w1 - w0]; %#ok<AGROW>
                    end
                end
                if ~isempty(Cp) && ~isempty(Tp), rowres = max(rowres, norm(Cp - (Cp*pinv(Tp))*Tp, 'fro')/norm(Cp, 'fro')); end
            end
        end
    end
    m = size(Arows, 1); nse = noise_ns*randn(m, 1);
    AA = Arows'*Arows; lams = logspace(-1, 4, 50);
    orac = @(dd) min(arrayfun(@(L) sqrt(mean((1./(1/c0 + ((AA + L^2*GG) \ (Arows'*dd))*1e-6) - ct).^2)), lams));
    rmse_cons = orac(Arows*dsf(ct) + nse);          % data consistent with the cell model: noise + regularization only
    rmse_nf = orac(d);                              % fine-ray data without noise: model error + regularization only
    dn = d + nse; Ab = Arows'*dn; rr = zeros(size(lams));
    for q = 1:numel(lams), sq = (AA + lams(q)^2*GG) \ Ab; rr(q) = sqrt(mean((1./(1/c0 + sq*1e-6) - ct).^2)); end
    [~, qd] = min(rr); lam = lams(qd);              % oracle lambda (minimum RMSE against the truth)
    ds_hat = (AA + lam^2*GG) \ Ab; c_hat = 1./(1/c0 + ds_hat*1e-6);
    Rm = (AA + lam^2*GG) \ AA;
    inc = hypot(reshape(XX', [], 1) - 2, reshape(ZZ', [], 1) - 20) < 3;
    fat = reshape(ZZ', [], 1) < 5;
    ixc = find(abs(xc - 1.2) < 1e-9);
    fwa = zeros(1, 3); fwl = zeros(1, 3); zq = [9 19 31];
    for t = 1:3
        iz = find(abs(zc - zq(t)) < 1.3, 1); k = (iz - 1)*Nx + ixc;
        row = reshape(Rm(k, :), Nx, Nz)';           % Nz x Nx kernel
        fwa(t) = fw(max(row(:, ixc), 0)', h); fwl(t) = fw(max(row(iz, :), 0), h);
    end
    fprintf('  %-8s %5d   %6.2f   %7.1f    %7.1f      %7.1f                %7.1f            %4.1f/%3.1f  %4.1f/%3.1f  %4.1f/%3.1f\n', designs{dsg}, m, lam, trace(Rm), ...
        sqrt(mean((c_hat - ct).^2)), mean(c_hat(inc)), mean(c_hat(fat)), fwa(1), fwl(1), fwa(2), fwl(2), fwa(3), fwl(3));
    fprintf('            RMSE with cell-consistent data + noise: %.1f m/s; fine-ray data without noise: %.1f m/s\n', rmse_cons, rmse_nf);
    sv = svd(Arows); res(dsg).sv = sv/sv(1); res(dsg).c = c_hat; res(dsg).R = Rm;
end
fprintf('Singular values above 1e-2 of the largest: tx-only %d, CMA %d (of %d unknowns)\n', sum(res(1).sv > 1e-2), sum(res(2).sv > 1e-2), np);
fprintf('CMA rows of each pixel lie in the span of its tx-only rows: max relative residual %.1e\n', rowres);

figure('Position', [30 30 1150 640]); colormap(jet);
subplot(2, 3, 1); imagesc(xc, zc, reshape(ct, Nx, Nz)'); axis image; caxis([1440 1590]); colorbar; title('truth (m/s)'); ylabel('z (mm)');
subplot(2, 3, 2); imagesc(xc, zc, reshape(res(1).c, Nx, Nz)'); axis image; caxis([1440 1590]); colorbar; title('tx-only, Tikhonov');
subplot(2, 3, 3); imagesc(xc, zc, reshape(res(2).c, Nx, Nz)'); axis image; caxis([1440 1590]); colorbar; title('CMA, Tikhonov');
k = (find(abs(zc - 19) < 1.3, 1) - 1)*Nx + find(abs(xc - 1.2) < 1e-9);
subplot(2, 3, 4); semilogy(res(1).sv, 'b', 'LineWidth', 1.4); hold on; semilogy(res(2).sv, 'r', 'LineWidth', 1.4); grid on;
xlabel('index'); ylabel('\sigma_i / \sigma_1'); legend('tx-only', 'CMA', 'Location', 'southwest'); title('singular spectrum of A');
subplot(2, 3, 5); imagesc(xc, zc, reshape(res(1).R(k, :), Nx, Nz)'); axis image; colorbar; title('R row at (1.2, 19.4) mm, tx-only'); xlabel('x (mm)');
subplot(2, 3, 6); imagesc(xc, zc, reshape(res(2).R(k, :), Nx, Nz)'); axis image; colorbar; title('R row at (1.2, 19.4) mm, CMA'); xlabel('x (mm)');
