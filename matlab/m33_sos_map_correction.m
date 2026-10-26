%% m33_sos_map_correction.m
% Correction from a sound-speed map (the Module 32 route): delays computed
% through the map versus one measured law, and the price of straight rays,
% a blurred map and a biased layer speed.
% Purpose : 128 elements (pitch 0.3 mm), 5 MHz, beamforming speed c0 = 1540
%           m/s. A fat layer (1450 m/s) runs from the skin to a wavy interface
%           h(x) = 12 + 1.5 sin(2 pi x / 15 mm) mm over 1540 m/s tissue.
%           True one-way times from element x_u to a point (x_p, z): Fermat's
%           principle with one refracting interface (1D search over the
%           crossing point). The aberration law is the excess over c0:
%             tau_u(x_p) = T_true(u, p) - r_up / c0.
%           Corrections compared at points x_p = -12 ... 12 mm, z = 40 mm:
%             none; one law measured at x_p = 0 and used everywhere;
%             straight-ray integration of the true map;
%             straight rays through a blurred map (Gaussian, 3 mm FWHM);
%             straight rays with the fat speed 10 m/s too high.
%           Residual = law - correction, mean and ramp removed; Strehl at 5 MHz.
% Change  : c_fat, amp, period, blur_mm, dc_bias, z.
% Expect  : the single law fails a few mm off the guide point; the
%           straight-ray map leaves a few ns (refraction ignored), the blurred
%           map more, the biased speed a smooth curvature error.
% Runtime : about 5 s in Octave.
clear; close all;

c0 = 1540; c_fat = 1450; f0 = 5e6; N = 128; p = 0.3e-3; xe = ((0:N-1) - (N-1)/2)*p;
amp = 1.5e-3; per = 15e-3; h = @(x) 12e-3 + amp*sin(2*pi*x/per);
z = 40e-3; xp_list = (-12:1:12)*1e-3; blur_mm = 3; dc_bias = 10;
xi = (-30e-3:0.02e-3:30e-3);                                                    % candidate crossing points
detr = @(v) v(:) - [ones(numel(v), 1) xe(:)]*([ones(numel(v), 1) xe(:)]\v(:));
Ss = @(v) abs(mean(exp(1j*2*pi*f0*v)))^2;

% gridded map for straight-ray integration
dg = 0.1e-3; xg = -30e-3:dg:30e-3; zg = 0:dg:z; [XG, ZG] = meshgrid(xg, zg);
s_true = 1/c0*ones(size(XG)); s_true(ZG < h(XG)) = 1/c_fat;
sig = blur_mm*1e-3/2.355; kx = -4*sig:dg:4*sig; gk = exp(-kx.^2/(2*sig^2)); gk = gk/sum(gk);
s_blur = conv2(conv2(s_true, gk, 'same'), gk', 'same');                         % blurred map (tomography-like resolution)
s_blur(:, [1:40, end-39:end]) = s_true(:, [1:40, end-39:end]);                  % avoid zero-padding at the lateral edges
s_blur(1:40, :) = s_true(1:40, :); s_blur(end-39:end, :) = s_true(end-39:end, :);
s_bias = 1/c0*ones(size(XG)); s_bias(ZG < h(XG)) = 1/(c_fat + dc_bias);
tt = linspace(0, 1, 801)';                                                      % samples along each straight ray
ray = @(S, xp) sum(interp2(xg, zg, S, xe + (xp - xe).*tt, z*tt*ones(1, N), 'linear'))/numel(tt).*hypot(xp - xe, z);   % straight-ray times, all elements

law = @(xp) arrayfun(@(xu) min(hypot(xi - xu, h(xi))/c_fat + hypot(xp - xi, z - h(xi))/c0), xe) - hypot(xp - xe, z)/c0;   % Fermat
L0 = law(0);
res = zeros(numel(xp_list), 5);
for q = 1:numel(xp_list)
    xp = xp_list(q); Lt = law(xp);
    ts = ray(s_true, xp) - hypot(xp - xe, z)/c0;
    tb = ray(s_blur, xp) - hypot(xp - xe, z)/c0;
    tbi = ray(s_bias, xp) - hypot(xp - xe, z)/c0;
    R = [detr(Lt), detr(Lt - L0), detr(Lt - ts), detr(Lt - tb), detr(Lt - tbi)];
    res(q, :) = std(R);
    for m = 1:5, St(q, m) = Ss(R(:, m)); end
end
names = {'none', 'one law (x=0)', 'straight-ray map', 'blurred map', 'biased fat speed'};
fprintf('Fat 1450 m/s to h(x) = 12 + 1.5 sin(2 pi x/15) mm, z = 40 mm, 5 MHz, 128 elements\n');
fprintf('Residual rms (ns, mean and ramp removed) / Strehl\n   x_p   ');
fprintf('%-20s', names{:}); fprintf('\n');
for q = 1:3:numel(xp_list)
    fprintf('%5.0f mm ', xp_list(q)*1e3);
    for m = 1:5, fprintf('%6.1f / %.2f        ', res(q, m)*1e9, St(q, m)); end
    fprintf('\n');
end
fprintf('\nMean over x_p = -12..12 mm:  '); for m = 1:5, fprintf('%s %.1f ns (S %.2f);  ', names{m}, mean(res(:, m))*1e9, mean(St(:, m))); end; fprintf('\n');
ii = find(St(:, 2) >= 0.8); fprintf('One law keeps S >= 0.8 for |x_p| <= %.0f mm\n', max(abs(xp_list(ii)))*1e3);

figure('Position', [50 50 900 330]);
subplot(1, 2, 1); imagesc(xg*1e3, zg*1e3, 1./s_true); axis image; colorbar; hold on; plot([xe(1) 0 xe(end)]*1e3, [0 z 0]*1e3, 'w-'); xlabel('x (mm)'); ylabel('z (mm)'); title('sound-speed map (m/s)');
subplot(1, 2, 2); plot(xp_list*1e3, St(:, 1), 'k', xp_list*1e3, St(:, 2), 'r', xp_list*1e3, St(:, 3), 'b', xp_list*1e3, St(:, 4), 'b--', xp_list*1e3, St(:, 5), 'g');
xlabel('lateral position of the point (mm)'); ylabel('Strehl after correction'); legend(names, 'Location', 'south'); ylim([0 1.05]);
