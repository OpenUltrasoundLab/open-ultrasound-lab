%% m26_compounding_speckle.m  --  incoherent vs coherent compounding of speckle looks
% Purpose : N looks of the same scatterers from sub-apertures of width w whose centres are
%           translated by d (all focused at the same point). In k-space each look keeps the
%           same axial band and a lateral coarray triangle shifted by 2d (in aperture units).
%           Field correlation between looks = overlap of the two triangles:
%             rho_F(s) = 1 - 1.5 s^2 + 0.75 s^3 (s <= 1),  (2 - s)^3 / 4 (1 < s < 2),  s = 2 d / w
%           Intensity correlation rho_I = rho_F^2. Summing intensities of N looks gives
%             C_I = sqrt(sum_ij rho_I(i,j)) / N,   N_eff = N^2 / sum_ij rho_I(i,j)
%           Coherent summation (IQ) widens the k-space support: finer speckle, but it is
%           still circular Gaussian, so envelope SNR0 stays 1.91 for any N.
% Change  : Nlist, dwlist (d/w), grid.
% Expect  : d/w >= 0.5 makes looks almost independent (rho_I <= 0.06): C_I ~ 1/sqrt(N);
%           coherent SNR0 = 1.91 +- 0.03 for every case; coherent lateral speckle shrinks.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);
Nz = 128;  Nx = 512;                           % grid (pixels); lateral unit = 1 pixel
ux = [0:Nx/2-1, -Nx/2:-1] / Nx;                % lateral frequency (cycles/pixel)
uz = [0:Nz/2-1, -Nz/2:-1]' / Nz;
U  = 1/32;                                     % coarray half-width of one look (cycles/pixel)
Hz = exp(-uz.^2 / (2*0.08^2));                 % axial band (Gaussian)
tri = @(u, c0) max(0, 1 - abs(u - c0)/U);
rhoF = @(s) (abs(s) <= 1).*(1 - 1.5*s.^2 + 0.75*abs(s).^3) + (abs(s) > 1 & abs(s) < 2).*(2 - abs(s)).^3/4;
acfx = @(I) mean(real(ifft(abs(fft(bsxfun(@minus, I, mean(I, 2)), [], 2)).^2, [], 2)), 1);
fw = @(p) 2*((find(p/p(1) < 0.5, 1) - 2) + (p(find(p/p(1) < 0.5, 1) - 1)/p(1) - 0.5) / ...
          ((p(find(p/p(1) < 0.5, 1) - 1) - p(find(p/p(1) < 0.5, 1)))/p(1)));
W = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
S = fft2(W);
single_look = ifft2(S .* (Hz * tri(ux, 0)));
w1 = fw(acfx(abs(single_look).^2));

Nlist = [1 3 5 9];  dwlist = [0.1 0.25 0.5 1.0];
fprintf('single look: SNR0 %.3f, lateral speckle FWHM %.1f px\n', mean(abs(single_look(:)))/std(abs(single_look(:))), w1);
fprintf('  N  d/w  | incoh. intensity C_I  meas / theory  N_eff | incoh. env SNR0 | coherent SNR0  lat. FWHM (px)\n');
res = zeros(numel(Nlist), numel(dwlist), 3);
for a = 1:numel(Nlist)
    N = Nlist(a);
    for b = 1:numel(dwlist)
        s = 2*dwlist(b);  ctr = ((0:N-1) - (N-1)/2) * s * U;
        Isum = 0;  Asum = 0;  Zsum = 0;
        for k = 1:N
            zk = ifft2(S .* (Hz * tri(ux, ctr(k))));
            Isum = Isum + abs(zk).^2;  Asum = Asum + abs(zk);  Zsum = Zsum + zk;
        end
        [ii, jj] = meshgrid(1:N);
        rhoI = rhoF((ii - jj)*s).^2;
        Cth = sqrt(sum(rhoI(:)))/N;  Neff = N^2/sum(rhoI(:));
        Cm = std(Isum(:))/mean(Isum(:));
        snrA = mean(Asum(:))/std(Asum(:));
        snrC = mean(abs(Zsum(:)))/std(abs(Zsum(:)));
        wc = fw(acfx(abs(Zsum).^2));
        res(a, b, :) = [Cm Cth snrC];
        fprintf('%3d  %4.2f |      %5.3f / %5.3f    %5.2f |     %5.3f       |    %5.3f       %5.1f\n', ...
            N, dwlist(b), Cm, Cth, Neff, snrA, snrC, wc);
        if N == 5 && dwlist(b) == 0.5, Iinc = Asum; Icoh = abs(Zsum); end
    end
end
fprintf('limit for independent looks: C_I = 1/sqrt(N); envelope SNR0 -> 1.91*sqrt(N)\n');

%% Figure
figure('Position', [60 60 1150 420]);
show = @(A) 20*log10(A/mean(A(:)));
xr = 1:256;
subplot(2,3,1); imagesc(show(abs(single_look(:, xr))), [-25 8]); colormap(gray); axis image; title('single look'); xlabel('lateral (px)'); ylabel('depth (px)');
subplot(2,3,2); imagesc(show(Iinc(:, xr)), [-25 8]); axis image; title('incoherent, N = 5, d/w = 0.5'); xlabel('lateral (px)');
subplot(2,3,3); imagesc(show(Icoh(:, xr)), [-25 8]); axis image; title('coherent, N = 5, d/w = 0.5'); xlabel('lateral (px)');
subplot(2,3,[4 5 6]);
cols = 'brkm';
for a = 1:numel(Nlist)
    plot(dwlist, squeeze(res(a,:,1)), [cols(a) 'o'], 'MarkerFaceColor', cols(a)); hold on;
    dg = linspace(0, 1, 101);  th = zeros(size(dg));
    for q = 1:numel(dg)
        [ii, jj] = meshgrid(1:Nlist(a));  r = rhoF((ii - jj)*2*dg(q)).^2;  th(q) = sqrt(sum(r(:)))/Nlist(a);
    end
    plot(dg, th, [cols(a) '-']);
end
xlabel('look translation d / sub-aperture width w'); ylabel('intensity speckle contrast C_I');
title('incoherent compounding, N = 1, 3, 5, 9: simulation (o) vs sqrt(sum rho_I)/N'); grid on; ylim([0 1.1]);
