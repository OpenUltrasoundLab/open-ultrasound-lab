%% m19_microbeamforming.m  --  sub-array (micro-) beamforming: subarray factor x coarse array factor
% Purpose : a 48-element line (one axis of a 48 x 48 matrix) at pitch d is grouped into
%           subarrays of m elements. Fine delays inside each subarray steer it to u0 = sin(th0)
%           (pre-steer); coarse delays between subarray centres steer to us = sin(ths). The
%           receive pattern factorizes exactly: P(u) = CAF(u - us) * SF(u - u0), where CAF has
%           grating lobes at us +- p lambda/(m d) and SF has nulls at u0 +- p lambda/(m d).
%           (1) checks the factorization against the direct element sum; (2) tabulates the
%           main-lobe loss SF(delta) and the highest coarse grating lobe versus the pre-steer
%           error delta = ths - th0 for m = 2, 3, 4; (3) the fine-delay span and the curvature
%           error inside a subarray; (4) 2D patterns of the 48 x 48 matrix (separable).
% Change  : m, d/lambda (dl), delta values, N, R.
% Expect  : factorization error ~1e-15. d = lambda/2, 3 x 3: delta = 5 deg: -0.22 dB, lobe -15.2 dB;
%           10 deg: -0.88 dB, -8.3 dB; 20 deg: -3.73 dB, lobe above the main lobe (+0.4 dB).
%           2 x 2 at 10 deg: -0.33 dB, -11.1 dB. 4 x 4 at 10 deg: -1.69 dB, -5.2 dB.
%           Channels of a 48 x 48 matrix: 2304 -> 256 (3 x 3). Fine-delay span (m-1) d sin(th0)/c
%           for 3 x lambda/2 at 37 deg, 5 MHz: 120 ns. Curvature phase pi (m d/2)^2/(lambda R) at
%           R = 20 mm: 0.027 rad.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

c = 1540; f0 = 5e6; lam = c/f0;                     % m/s, Hz, m
N = 48; dl = 0.5; d = dl*lam;                        % elements, pitch in lambda, pitch (m)
u = linspace(-1, 1, 4001);
SF = @(u, m) abs(sin(m*pi*dl*u) ./ (m*sin(pi*dl*u) + (abs(sin(pi*dl*u)) < 1e-12))) .* (abs(sin(pi*dl*u)) >= 1e-12) + (abs(sin(pi*dl*u)) < 1e-12);

%% (1) Direct element sum versus CAF x SF
m = 3; K = N/m; u0 = sin(0*pi/180); us = sin(10*pi/180);
n = 0:N-1; k = floor(n/m); xi = mod(n, m) - (m-1)/2; Xk = (k + 0.5)*m - N/2;   % positions in pitch units
w = exp(-1i*2*pi*dl*(xi*u0 + Xk*us));                                        % fine + coarse steering
Pdir = abs(exp(1i*2*pi*dl*u(:)*(Xk + xi)) * w(:)).' / N;
CAF  = abs(exp(1i*2*pi*dl*(u(:) - us)*(((0:K-1) + 0.5)*m - N/2)) * ones(K, 1)).' / K;
Pfac = CAF .* SF(u - u0, m);
fprintf('(1) max |direct - CAF*SF| = %.1e  (m = %d, delta = 10 deg)\n', max(abs(Pdir - Pfac)), m);

%% (2) Loss and grating lobe versus pre-steer error (th0 = 0)
dels = [0 5 10 15 20];
fprintf('\n(2) d = %.2f lambda: main-lobe loss / highest coarse grating lobe re main (dB)\n   delta', dl);
fprintf('       m = %d      ', 2:4); fprintf('\n');
for dd = dels
    fprintf('   %4.0f   ', dd);
    for m = 2:4
        ust = sin(dd*pi/180); loss = 20*log10(SF(ust, m)); gl = -Inf;
        for p = [-2*m:-1, 1:2*m]
            if mod(p, m) == 0, continue; end
            ug = ust + p/(m*dl);
            if abs(ug) <= 1, gl = max(gl, 20*log10(SF(ug, m) + 1e-15) - loss); end
        end
        if gl < -100, fprintf('%6.2f / none    ', loss); else, fprintf('%6.2f / %+5.1f   ', loss, gl); end
    end
    fprintf('\n');
end
fprintf('   channels for 48 x 48: %d -> %d (2x2), %d (3x3), %d (4x4)\n', N^2, (N/2)^2, (N/3)^2, (N/4)^2);

%% (3) Fine-delay span and curvature error
th0 = 37*pi/180; R = 20e-3;
fprintf('\n(3) fine-delay span (m-1) d sin(th0)/c, 3 x lambda/2 at 37 deg: %.0f ns\n', 2*d*sin(th0)/c*1e9);
for mm = [3 4 8]
    fprintf('    curvature phase pi (m d/2)^2/(lambda R), m = %d: R = 5 mm %.3f rad, R = 20 mm %.3f rad\n', ...
        mm, pi*(mm*d/2)^2/(lam*5e-3), pi*(mm*d/2)^2/(lam*R));
end

%% (4) Figures: 1D cuts and 2D patterns for delta = 0 and 10 deg (m = 3)
m = 3; K = N/m; k = floor(n/m); xi = mod(n, m) - (m-1)/2; Xk = (k + 0.5)*m - N/2;
Pof = @(u0, us) abs(exp(1i*2*pi*dl*u(:)*(Xk + xi)) * exp(-1i*2*pi*dl*(xi*u0 + Xk*us)).').' / N;
P0 = Pof(0, 0); P10 = Pof(0, sin(10*pi/180));
figure('Name', 'Micro-beamforming cuts');
subplot(2, 1, 1); plot(u, 20*log10(P0 + 1e-6), 'b', u, 20*log10(SF(u, m) + 1e-6), 'r--'); ylim([-50 0]); grid on;
title('3 x 3 at \lambda/2, beam = pre-steer (\delta = 0)'); ylabel('dB');
subplot(2, 1, 2); plot(u, 20*log10(P10 + 1e-6), 'b', u, 20*log10(SF(u, m) + 1e-6), 'r--'); ylim([-50 0]); grid on;
title('\delta = 10 deg: subarray factor no longer nulls the coarse lobes'); xlabel('u = sin \theta'); ylabel('dB');
[UU, VV] = meshgrid(u(1:20:end), u(1:20:end));
vis = (UU.^2 + VV.^2) <= 1;
img = @(Pu) max(20*log10(Pu(1:20:end).' * P0(1:20:end) + 1e-9), -40) .* vis - 45*(~vis);
figure('Name', 'Micro-beamforming 2D');
subplot(1, 2, 1); imagesc(u(1:20:end), u(1:20:end), img(P0).');  axis image; caxis([-40 0]); xlabel('u'); ylabel('v'); title('\delta = 0');
subplot(1, 2, 2); imagesc(u(1:20:end), u(1:20:end), img(P10).'); axis image; caxis([-40 0]); xlabel('u'); title('\delta = 10 deg in u'); colorbar;
