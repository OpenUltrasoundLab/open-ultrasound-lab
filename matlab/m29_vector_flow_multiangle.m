%% m29_vector_flow_multiangle.m
% Purpose : 2D velocity vectors from multi-angle plane-wave Doppler.
%           Each transmit/receive angle pair (alpha_k, beta_k) measures the
%           projection f_k = -(f0/c) [ (sin a + sin b) vx + (cos a + cos b) vz ]
%           (x lateral, z depth downward). K >= 2 pairs give a least-squares
%           vector. Part 1: Monte Carlo (Kasai on simulated ensembles) versus
%           the covariance sigma_f^2 (A'A)^-1. Part 2: the cross-beam cot(gamma)
%           law for the lateral/axial precision ratio. Part 3: the lateral
%           wavelength of transverse oscillation, lambda_x = 2 lambda z / d,
%           measured from the receive field of a two-peak apodization.
% Change  : gamma (angle span), K, mode (rx = tx or rx = 0), snrdB, N.
% Expect  : precision ratio sigma_vx/sigma_vz = cot(gamma) for 2 symmetric
%           pairs with rx = tx: 5.67 at 10 deg, 2.75 at 20 deg. With rx fixed
%           at 0 deg the ratio is cot(gamma/2), about twice as large.
%           lambda_x: z = 30 mm, d = 10 mm, f0 = 5 MHz -> 1.85 mm.
rng(1);

%% Parameters
c = 1540; f0 = 5e6; PRFtx = 10e3;   % transmit rate (Hz); per-angle PRF = PRFtx/K
N = 32;  snrdB = 10;  nT = 2000;    % ensemble per angle, per-sample SNR, trials
v = [0.20, -0.10];                  % true (vx, vz) in m/s; vz < 0 moves toward the probe
Arow = @(al, be) -(f0/c) * [sin(al) + sin(be), cos(al) + cos(be)];

%% Part 1 and 2: Monte Carlo vs covariance
fprintf('True v = (%.2f, %.2f) m/s, N = %d per angle, SNR = %d dB, per-angle PRF = PRFtx/K\n', v(1), v(2), N, snrdB);
fprintf(' mode     K  gamma  sd_vx MC  sd_vx th  sd_vz MC  sd_vz th  ratio  cot(g) cot(g/2)  v_Nyq,axial(m/s)\n');
modes = {'rx=tx', 'rx=0'};
for im = 1:2
    for K = [2 5]
        for g = [5 10 20 30]
            gam = g*pi/180;  al = linspace(-gam, gam, K)';  be = al * (im == 1);
            A = zeros(K, 2);  for k = 1:K, A(k, :) = Arow(al(k), be(k)); end
            PRF = PRFtx/K;  f = A*v';  m = (0:N-1)';  est = zeros(nT, 2);  fh = zeros(K, nT);
            for k = 1:K
                z = exp(2j*pi*f(k)*m/PRF) * exp(2j*pi*rand(1, nT)) + (randn(N, nT) + 1j*randn(N, nT))/sqrt(2*10^(snrdB/10));
                fh(k, :) = PRF/(2*pi) * angle(sum(z(2:end, :).*conj(z(1:end-1, :)), 1));
            end
            est = (A \ fh)';                                 % least squares, one vector per trial
            sf = std(fh(1, :) - f(1));                      % per-angle Doppler std (Hz)
            Cv = sf^2 * inv(A'*A);
            fprintf(' %-6s  %d   %2d    %.4f    %.4f    %.4f    %.4f   %5.2f  %5.2f  %5.2f    %.3f\n', modes{im}, K, g, std(est(:,1)), sqrt(Cv(1,1)), ...
                std(est(:,2)), sqrt(Cv(2,2)), std(est(:,1))/std(est(:,2)), cot(gam), cot(gam/2), c*PRF/(4*f0));
            if im == 1 && K == 2 && g == 10, cloud = est; end
        end
    end
end
fprintf('Mean estimate (rx=tx, K = 2, 10 deg): (%.3f, %.3f) m/s\n', mean(cloud(:,1)), mean(cloud(:,2)));

%% Part 3: transverse oscillation lateral wavelength from a two-peak receive apodization
lam = c/f0;  d = 10e-3;  w = 2.5e-3;           % peak separation and Gaussian peak width (m)
xi = linspace(-20e-3, 20e-3, 4096);  dxi = xi(2) - xi(1);
apod = exp(-((xi - d/2)/w).^2) + exp(-((xi + d/2)/w).^2);
fprintf('Transverse oscillation, f0 = %.0f MHz, d = %.0f mm:\n   z (mm)  lambda_x measured (mm)  2*lambda*z/d (mm)\n', f0/1e6, d*1e3);
for zf = [20 30 40]*1e-3
    x = linspace(-6e-3, 6e-3, 1201);
    h = exp(-2j*pi*x'*xi/(lam*zf)) * apod.' * dxi;   % Fraunhofer receive field at depth zf (focused)
    hr = real(h);                                    % field is real (symmetric apodization)
    big = abs(x(1:end-1)') < 3*lam*zf/(pi*w);        % inside 3x the Gaussian envelope radius (zeros of h are the crossings, so do not mask on |h|)
    zc = find(big & (hr(1:end-1) > 0 & hr(2:end) <= 0 | hr(1:end-1) < 0 & hr(2:end) >= 0));
    pz = mean(diff(x(zc))) * 2;                     % period from zero crossings
    fprintf('   %4.0f        %6.3f                %6.3f\n', zf*1e3, pz*1e3, 2*lam*zf/d*1e3);
    if abs(zf - 30e-3) < 1e-6, xTO = x; hTO = h; end
end

%% Figures
figure('Position', [50 50 1100 360]);
subplot(1,3,1); plot(cloud(:,1), cloud(:,2), '.', 'MarkerSize', 3); hold on; plot(v(1), v(2), 'k+', 'MarkerSize', 14, 'LineWidth', 2);
axis equal; grid on; xlabel('v_x (m/s)'); ylabel('v_z (m/s)'); title('LS estimates, K = 2, \pm10 deg, rx = tx');
gg = linspace(2, 40, 100)*pi/180;
subplot(1,3,2); semilogy(gg*180/pi, cot(gg), 'LineWidth', 1.5); hold on; semilogy(gg*180/pi, cot(gg/2), '--', 'LineWidth', 1.5);
grid on; xlabel('\gamma (deg)'); ylabel('\sigma_{vx} / \sigma_{vz}'); legend('rx = tx: cot \gamma', 'rx = 0: cot(\gamma/2)'); title('Lateral precision penalty');
subplot(1,3,3); plot(xTO*1e3, real(hTO)/max(abs(hTO)), xTO*1e3, abs(hTO)/max(abs(hTO)), '--', 'LineWidth', 1.3);
grid on; xlabel('lateral x (mm)'); ylabel('receive field (norm.)'); title('TO field, z = 30 mm, d = 10 mm');
