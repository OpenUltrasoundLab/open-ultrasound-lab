%% m34_persistence.m  --  frame persistence: time constant, noise reduction, motion lag, frame rate
% Purpose : First-order recursive (IIR) persistence  y_n = (1 - alpha) y_{n-1} + alpha x_n
%           at frame period T_F. Time constant tau = -T_F / ln(1 - alpha); mean lag
%           (1 - alpha)/alpha frames; for independent frames the noise variance falls by
%           alpha/(2 - alpha), i.e. N_eff = (2 - alpha)/alpha frames. For frames whose noise has
%           frame-to-frame correlation rho^|k| (AR(1)), with beta = 1 - alpha:
%             var(y)/var(x) = alpha/(2 - alpha) * (1 + beta rho)/(1 - beta rho).
%           Static speckle (rho = 1) is not reduced at all. A fixed alpha means a time constant
%           that changes with frame rate; a fixed tau needs alpha = 1 - exp(-T_F/tau).
%           Last part: a point target moving laterally, imaged at f_F, displayed with persistence:
%           the displayed peak drops and its centroid lags by v (1 - alpha)/alpha T_F.
% Change  : frame rates fF, alpha, tau, rho list, target speed v, beam width.
% Expect  : alpha = 0.5 gives tau = 72 ms at 20 Hz and 14 ms at 100 Hz; tau = 100 ms needs
%           alpha = 0.39 at 20 Hz and 0.095 at 100 Hz (N_eff 4.1 and 20). Monte Carlo variance
%           ratios match the AR(1) formula to about 1 %.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(5);

%% 1. alpha, tau and frame rate
fF = [20 30 60 100 200];
fprintf('fixed alpha = 0.5:   ');  fprintf('%5.0f Hz: tau %5.1f ms | ', [fF; -1e3./fF/log(0.5)]);  fprintf('\n');
tau = 0.1;  al = 1 - exp(-1./(fF*tau));
fprintf('fixed tau = 100 ms:  ');  fprintf('%5.0f Hz: alpha %.3f, N_eff %4.1f | ', [fF; al; (2 - al)./al]);  fprintf('\n');
fprintf('mean lag (1-alpha)/alpha * T_F at tau = 100 ms: ');  fprintf('%.0f ', 1e3*(1 - al)./al./fF);  fprintf('ms\n');

%% 2. noise reduction vs frame-to-frame correlation (Monte Carlo)
alpha = [0.5 0.25 0.1];  rho = [0 0.5 0.9 1];  Nf = 4000;  Np = 400;
fprintf('\n alpha  rho  | var ratio MC   AR(1) theory\n');
for a = alpha
    for r = rho
        x = zeros(Np, Nf);  x(:, 1) = randn(Np, 1);
        for n = 2:Nf, x(:, n) = r*x(:, n-1) + sqrt(1 - r^2)*randn(Np, 1); end
        y = zeros(Np, Nf);  y(:, 1) = x(:, 1);
        for n = 2:Nf, y(:, n) = (1 - a)*y(:, n-1) + a*x(:, n); end
        b = 1 - a;  th = a/(2 - a)*(1 + b*r)/(1 - b*r);
        v = var(reshape(y(:, 501:end), [], 1)) / var(reshape(x(:, 501:end), [], 1));
        fprintf(' %4.2f  %4.2f |   %6.3f       %6.3f\n', a, r, v, th);
    end
end

%% 3. moving point target: displayed profile with persistence
wb = 1.0e-3;                                    % lateral beam -6 dB width (m)
s = wb/(2*sqrt(2*log(2)));  xg = (-45:0.02:10)*1e-3;  v = 20e-3;   % 20 mm/s lateral motion
fprintf('\n target 20 mm/s, beam 1.0 mm: frame rate, alpha -> displayed peak, lag (mm), width (mm)\n');
for cfg = [20 0.5; 20 0.2; 100 0.5; 100 1 - exp(-1/(100*0.1))]'
    fr = cfg(1);  a = cfg(2);  T = 1/fr;  n = round(2/T);  y = zeros(size(xg));
    for k = 1:n
        xt = v*(k*T - n*T);                     % target position, ends at x = 0
        y = (1 - a)*y + a*exp(-(xg - xt).^2/(2*s^2));
    end
    [pk, ip] = max(y);  cen = sum(xg.*y)/sum(y);  hm = xg(y >= pk/2);
    fprintf('  %3.0f Hz, alpha %.3f (tau %5.1f ms): peak %.2f, centroid lag %.2f mm (theory %.2f), FWHM %.2f mm\n', ...
            fr, a, -1e3*T/log(1 - a), pk, -cen*1e3, v*(1 - a)/a*T*1e3, (hm(end) - hm(1))*1e3);
    if fr == 20 && a == 0.2, y20 = y; end
end

%% 4. frame rate limits from depth and line count
c = 1540;  zmax = [0.08 0.15 0.20];  Nl = [128 256];
fprintf('\nmax frame rate c/(2 z_max N_lines) (one transmit per line):\n');
for z = zmax, fprintf('  z_max %3.0f mm: ', z*1e3); fprintf('%4d lines %5.1f Hz  ', [Nl; c./(2*z*Nl)]); fprintf('\n'); end

figure('Position', [60 60 900 340]);
plot(xg*1e3, y20, 'r', xg*1e3, exp(-xg.^2/(2*s^2)), 'k--', 'LineWidth', 1.5); grid on;
xlabel('lateral position (mm)'); ylabel('displayed amplitude');
legend('20 Hz, alpha 0.2 (tau 224 ms)', 'target now (no persistence)', 'Location', 'northwest');
title('persistence smears a moving target into a trail behind it');
