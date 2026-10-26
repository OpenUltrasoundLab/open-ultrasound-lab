%% m32_cute_phase_shift.m
% CUTE-style echo phase shifts between steered acquisitions (ray model) and
% the phase jitter that limits them (speckle Monte Carlo).
% Purpose : (A) Straight-ray, first-order model of the beamforming delay
%           error of every pixel r for a steering angle th:
%              L_th(r) = sec(th) * int_0^z ds(x - (z - z') tan th, z') dz',
%           ds = 1/c - 1/c0 (c0 = 1540 m/s assumed by the beamformer).
%           Two data types are mapped for a +40 m/s inclusion under a fat layer:
%              transmit-only pair   d = L_th - L_0             (receive fixed)
%              common mid-angle     d = L_th + L_-th - 2 L_0   (tx th, rx -th)
%           Echo-position shifts and refraction are neglected.
%           (B) The phase of the complex correlation of two speckle kernels of
%           N independent samples and correlation rho has std
%              sigma_phi ~ sqrt(1 - rho^2) / (rho sqrt(2N))   (large N),
%           checked by Monte Carlo; sigma_tau = sigma_phi / (2 pi f0).
% Change  : th_deg, dc_inc, r_inc, the fat layer, rho_list, N_list.
% Expect  : at 15 deg the shift below the inclusion centre (0 deg path through
%           the full chord, 15 deg path missing it) reaches the chord delay
%           (164 ns) plus 5.7 ns from the fat layer, ~169 ns, several phase wraps at 5 MHz; a 2 deg step
%           stays below 100 ns (no wrap), zero above the inclusion; the layer adds only (sec th - 1) of its
%           vertical delay; the jitter formula holds within ~10 % for N = 64
%           and fails at low rho and small N (phase wrapping, bias).
% Runtime : about 2 s in Octave.
clear; close all; rng(1);

c0 = 1540; f0 = 5e6;
dx = 0.1; x = -19.2:dx:19.2; z = 0:dx:40;    % grid (mm)
[X, Z] = meshgrid(x, z);
dc_inc = 40; r_inc = 5; zi = 15;             % inclusion: +40 m/s, radius 5 mm, centre depth 15 mm
c = c0 + dc_inc./(1 + exp((hypot(X, Z - zi) - r_inc)/0.3));   % smooth edge, 0.3 mm
d_fat = 4; c(Z < d_fat) = 1450;              % 4 mm fat layer at the top
ds = (1./c - 1/c0)*1e6;                      % slowness deviation (ns/mm)

th_deg = 15; ths = [-th_deg 0 th_deg 2]*pi/180; L = cell(1, 4);
for a = 1:4
    th = ths(a); La = zeros(size(ds));
    La(1, :) = ds(1, :)*dx/cos(th)/2;        % half cell at the surface
    for k = 2:numel(z)                       % march down: shift the previous row by dz tan(th)
        prev = interp1(x, La(k-1, :), x - dx*tan(th), 'linear', 'extrap');
        La(k, :) = prev + 0.5*(ds(k, :) + ds(k-1, :))*dx/cos(th);
    end
    L{a} = La;                               % ns, delay error of each pixel for angle th
end
d_tx  = L{3} - L{2};                         % transmit-only pair (th vs 0)
d_cma = L{3} + L{1} - 2*L{2};                % common mid-angle pair (th,-th) vs (0,0)
in_ap = abs(X) < 12;                         % keep away from the lateral edges
fprintf('Inclusion %+d m/s: ds = %.2f ns/mm, vertical delay through its centre %.1f ns\n', dc_inc, ...
    (1/(c0 + dc_inc) - 1/c0)*1e6, (1/(c0 + dc_inc) - 1/c0)*1e6*2*r_inc);
fprintf('Peak |d| at th = %d deg: transmit-only %.1f ns (%.2f rad at %.0f MHz), common mid-angle %.1f ns (%.2f rad)\n', th_deg, ...
    max(abs(d_tx(in_ap))), 2*pi*f0*max(abs(d_tx(in_ap)))*1e-9, f0/1e6, max(abs(d_cma(in_ap))), 2*pi*f0*max(abs(d_cma(in_ap)))*1e-9);
d2 = L{4} - L{2};
fprintf('Adjacent-angle step of 2 deg (tx-only): peak |d| %.1f ns (%.2f rad): no wrapping; wrapping starts at |d| = 1/(2 f0) = %.0f ns\n', ...
    max(abs(d2(in_ap))), 2*pi*f0*max(abs(d2(in_ap)))*1e-9, 1e9/(2*f0));
k30 = find(abs(z - 30) < 1e-9); i0 = find(abs(x + 15) < 1e-9);   % a column far from the inclusion
S = (1/1450 - 1/c0)*1e6*d_fat;
fprintf('Fat layer only (x = -15 mm, z = 30 mm): vertical %.1f ns, tx-only %.2f ns (pred (sec-1)S = %.2f), CMA %.2f ns (pred 2(sec-1)S = %.2f)\n', ...
    L{2}(k30, i0), d_tx(k30, i0), (sec(ths(3)) - 1)*S, d_cma(k30, i0), 2*(sec(ths(3)) - 1)*S);

% (B) phase jitter of the complex correlation of two speckle kernels
rho_list = [0.5 0.7 0.9 0.95 0.99]; N_list = [16 64]; R = 4000;
fprintf('\n  rho     N   sigma_phi MC (rad)   formula   sigma_tau at 5 MHz (ns)\n');
SPH = zeros(numel(N_list), numel(rho_list));
for iN = 1:numel(N_list)
    N = N_list(iN);
    for ir = 1:numel(rho_list)
        rho = rho_list(ir); phi0 = 0.4;
        a1 = (randn(N, R) + 1j*randn(N, R))/sqrt(2);
        w  = (randn(N, R) + 1j*randn(N, R))/sqrt(2);
        a2 = (rho*a1 + sqrt(1 - rho^2)*w)*exp(1j*phi0);
        ph = angle(sum(a2.*conj(a1), 1)) - phi0;
        SPH(iN, ir) = std(ph);
        fo = sqrt(1 - rho^2)/(rho*sqrt(2*N));
        fprintf('  %.2f  %3d   %10.4f        %8.4f   %8.2f\n', rho, N, SPH(iN, ir), fo, SPH(iN, ir)/(2*pi*f0)*1e9);
    end
end
th = 20*pi/180; dzl = 5; sig_tau = sqrt(1 - 0.95^2)/(0.95*sqrt(128))/(2*pi*f0)*1e9;
% layer delay S = datum below minus datum above the layer: two independent data, sqrt(2) sigma_tau
sig_s = sqrt(2)*sig_tau/(2*(sec(th) - 1)*dzl);
fprintf('\nLayer slowness from CMA data below minus above a %d mm layer at +-20 deg, sigma_tau = %.2f ns per datum (rho 0.95, N 64): sigma_s = %.1f ns/mm = %.0f m/s\n', ...
    dzl, sig_tau, sig_s, sig_s*1e-6*c0^2);

figure('Position', [40 40 1150 380]); colormap(jet);
subplot(1, 3, 1); imagesc(x, z, c); axis image; colorbar; title('true c (m/s)'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(1, 3, 2); imagesc(x, z, d_tx); axis image; colorbar; caxis([-180 180]);
title(sprintf('tx-only: L_{%d} - L_0 (ns)', th_deg)); xlabel('x (mm)');
subplot(1, 3, 3); imagesc(x, z, d_cma); axis image; colorbar; caxis([-180 180]);
title(sprintf('CMA: L_{%d} + L_{-%d} - 2L_0 (ns)', th_deg, th_deg)); xlabel('x (mm)');
figure('Position', [60 60 560 380]);
semilogy(rho_list, SPH(1, :), 'bo', rho_list, SPH(2, :), 'ro', 'LineWidth', 1.4); hold on;
rr = linspace(0.45, 0.995, 200);
semilogy(rr, sqrt(1 - rr.^2)./(rr*sqrt(2*16)), 'b-', rr, sqrt(1 - rr.^2)./(rr*sqrt(2*64)), 'r-');
grid on; xlabel('correlation coefficient rho'); ylabel('phase std (rad)');
legend('MC N = 16', 'MC N = 64', 'formula N = 16', 'formula N = 64', 'Location', 'southwest');
title('Phase jitter of the complex correlation');
