%% m23_angle_step_rule.m  --  transmit replicas, the angle-step rule, and motion between angles
% Purpose : fast analytic plane-wave synthesis (no channel data): for a point target the
%           DAS value at pixel p for angle theta is
%             v(p, theta) = sum_e w_e(p) g(tau_tx(p,theta) + tau_rx(p,e) - tau_tx(t,theta) - tau_rx(t,e))
%           with the complex two-way pulse g(t) = exp(-t^2/(2 sig^2)) exp(j 2 pi f0 t).
%           (1) Replica position and level vs angle step (11 angles): the transmit factor
%               repeats every lambda/d(sin theta) laterally; we measure the broadband ghost
%               and compare it with the narrowband focal-plane prediction R(x_rep).
%           (2) Depth independence of the replica position (targets at 20, 30, 45 mm).
%           (3) Motion: the target moves axially between firings (31 angles, 1 deg step,
%               PRF 10 kHz). Loss at the true position vs the phase-ramp law
%               |sin(N phi/2)/(N sin(phi/2))|, phi = 2 pi f0 (1 + cos theta) delta / c per
%               firing; peak loss and lateral shift for sequential vs alternating order.
% Change  : steps, N, Fnum, FBW, PRF, velocities, the firing order.
% Expect  : the broadband ghost is a hump centered near lambda/d(sin theta) and smeared by
%           about +-FBW/2 of that distance, at the same x for every depth. Its level is far
%           below the narrowband R(x_rep), roughly R(x_rep) x 1.06 M/N: the virtual transmit
%           array has broadband grating lobes like a real array (-25 dB at a 15 deg step,
%           about -53 dB at 3 deg). Motion: at the true target position the loss follows the
%           ramp law (lambda/4 across the sequence = 25 mm/s here: -3.9 dB). In sequential
%           order the peak moves 0.12 mm laterally and drops only 2.2 dB (part of the loss
%           becomes a mislocation); in alternating order it stays put and loses 4.0 dB.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

c = 1.540; f0 = 5; lam = c/f0; FBW = 0.6;          % mm/us, MHz, mm
sig = 1/(2*pi*FBW*f0/(2*sqrt(2*log(2))));          % envelope std of the two-way pulse (us)
Nel = 128; pitch = 0.3; xe = ((0:Nel-1) - (Nel-1)/2)*pitch;
Fnum = 1.5;
g = @(t) exp(-t.^2/(2*sig^2)) .* exp(1i*2*pi*f0*t);
% geometry of a pixel set, computed once: receive times (npix x Nel) and normalized weights
geom = @(px, pz) struct('px', px, 'pz', pz, 'R', sqrt((px - xe).^2 + pz.^2)/c, ...
    'W', double(abs(xe - px) <= pz/(2*Fnum)) ./ sum(abs(xe - px) <= pz/(2*Fnum), 2));
% DAS value of one angle for a target at (xt, zt): sum_e W g(tx delay difference + rx delay difference)
das = @(G, xt, zt, th) sum(G.W .* g(((G.px - xt)*sin(th) + (G.pz - zt)*cos(th))/c ...
        + G.R - sqrt((xt - xe).^2 + zt^2)/c), 2);

%% (1) Replica vs angle step, 11 angles, target at (0, 30 mm)
zt = 30;
xl = (0:0.05:12).';  dzs = -0.6:0.1:0.6;                  % strip of pixels (x >= 0, symmetric)
[PX, DZ] = meshgrid(xl, dzs);  G1 = geom(PX(:), zt + DZ(:));
steps = [1 2 3 5 7.5 10 15];  N = 11;
% narrowband receive pattern at depth zt (paraxial), for the prediction
wz  = double(abs(xe) <= zt/(2*Fnum));
Rnb = @(x) abs(exp(-1i*2*pi*x(:)*xe/(lam*zt)) * wz.') / sum(wz);
fprintf('(1) 11 angles, target at z = %g mm, Rx F# %.1f, FBW %.0f %%\n', zt, Fnum, 100*FBW);
fprintf('%6s %5s %5s %8s %9s %8s %11s %11s\n', 'step', 'N', 'span', 'x_rep', 'hump at', 'ghost', 'narrowband', '+1.06M/N');
fprintf('%6s %5s %5s %8s %9s %8s %11s %11s\n', '(deg)', '', '(deg)', '(mm)', 'x (mm)', '(dB)', 'R(x_rep)', '(dB)');
M = 0.883/FBW;                                             % pulse length in cycles
prof = zeros(numel(xl), numel(steps));
for k = 1:numel(steps)
    if (N-1)/2*steps(k) > 25, N = 2*floor(25/steps(k)) + 1; end       % keep |theta| <= 25 deg
    th = ((0:N-1) - (N-1)/2) * steps(k) * pi/180;
    v = zeros(numel(G1.px), 1);
    for a = 1:N, v = v + das(G1, 0, zt, th(a)); end
    E = reshape(abs(v), size(PX));  E = max(E, [], 1).';  E = E/max(E);   % max over the strip
    prof(:,k) = 20*log10(E + 1e-9);
    xr = lam/mean(diff(sin(th)));
    win = xl > 0.6 & abs(xl - xr) <= 0.3 + 0.5*FBW*xr;       % the smear spans ~ +-FBW/2 of x_rep
    [gh, ig] = max(prof(win, k));  xw = xl(win);
    if xr > 11.5, gh = NaN; xm = NaN; else, xm = xw(ig); end
    fprintf('%6.1f %5d %5.0f %8.2f %9.2f %8.1f %11.1f %11.1f\n', steps(k), N, (N-1)*steps(k), xr, xm, gh, ...
        20*log10(Rnb(xr)), 20*log10(Rnb(xr)) + 20*log10(min(1, 1.064*M/N)));
    N = 11;
end

%% (2) The replica position does not depend on depth (3 deg step, 11 angles)
th = ((0:10) - 5) * 3 * pi/180;  xr = lam/mean(diff(sin(th)));
fprintf('\n(2) 3 deg step: predicted x_rep = %.2f mm at every depth\n', xr);
for ztk = [20 30 45]
    Gk = geom(PX(:), ztk + DZ(:));  v = zeros(numel(Gk.px), 1);
    for a = 1:numel(th), v = v + das(Gk, 0, ztk, th(a)); end
    E = max(reshape(abs(v), size(PX)), [], 1).';  E = E/max(E);
    win = xl > 3 & xl < 9;  xw = xl(win);  [gh, ig] = max(E(win));
    fprintf('    z = %2d mm: ghost hump peak at x = %.2f mm, level %.1f dB\n', ztk, xw(ig), 20*log10(gh));
end

%% (3) Motion between firings: 31 angles, 1 deg step, PRF 10 kHz, axial velocity v
th0 = (-15:15)*pi/180;  Na = numel(th0);  PRF = 10e3;          % firings per second
ordA = 1:Na;                                                   % sequential: -15, -14, ..., +15
ordB = [16 reshape([15:-1:1; 17:31], 1, [])];                 % alternating: 0, -1, +1, -2, +2, ...
xs = (-1:0.02:1).';  zs = zt + (-0.2:0.02:0.2);
[XS, ZS] = meshgrid(xs, zs);  G3 = geom(XS(:), ZS(:));  i0 = find(abs(XS(:)) < 1e-9 & abs(ZS(:) - zt) < 1e-9);
v0 = zeros(numel(G3.px), 1);
for a = 1:Na, v0 = v0 + das(G3, 0, zt, th0(a)); end
pk0 = max(abs(v0));
vel = [5 12.4 25 37];                                         % mm/s
fprintf('\n(3) motion, %d angles at PRF %.0f kHz (sequence %.1f ms); loss re static peak (dB)\n', Na, PRF/1e3, Na/PRF*1e3);
fprintf('%8s %11s %11s %11s %11s %13s %12s\n', 'v', 'disp', 'ramp law', 'seq: at', 'seq: peak', 'seq: shift', 'alt: peak');
fprintf('%8s %11s %11s %11s %11s %13s %12s\n', '(mm/s)', '(lambda)', '', 'target', 'anywhere', 'x (mm)', 'anywhere');
% disp = v x sequence duration = Na*dlt (the target moves dlt between firings)
res3 = zeros(numel(vel), 5);
for j = 1:numel(vel)
    dlt = vel(j)/PRF;                                         % mm per firing (v in mm/s)
    phi = 2*pi*f0*mean(1 + cos(th0))*dlt/c;                   % two-way phase step per firing
    law = 20*log10(abs(sin(Na*phi/2)/(Na*sin(phi/2))));
    vA = zeros(numel(G3.px), 1); vB = vA;
    for n = 1:Na                                              % firing n uses angle ord(n)
        zn = zt + (n - (Na+1)/2)*dlt;                         % target position at firing n
        vA = vA + das(G3, 0, zn, th0(ordA(n)));
        vB = vB + das(G3, 0, zn, th0(ordB(n)));
    end
    [pA, iA] = max(abs(vA));
    res3(j,:) = [law, 20*log10(abs(vA(i0))/pk0), 20*log10(pA/pk0), G3.px(iA), 20*log10(max(abs(vB))/pk0)];
    fprintf('%8.1f %11.3f %11.2f %11.2f %11.2f %13.2f %12.2f\n', vel(j), Na*dlt/lam, res3(j,:));
end
fprintf('    sequential order: the ramp is linear in sin(theta), so it shifts the transmit factor\n');
fprintf('    laterally by 2 dlt/d(sin theta) = %.2f mm at 25 mm/s against the fixed receive factor:\n', 2*(25/PRF)/mean(diff(sin(th0))));
fprintf('    the peak moves part of the way and loses less; the alternating order cannot shift.\n');

%% Figures
figure('Position', [60 60 1150 430]);
subplot(1,2,1); plot(xl, prof(:, [1 3 5 7]), 'LineWidth', 1.1); grid on; ylim([-80 3]); xlim([0 12]);
xlabel('lateral offset x (mm)'); ylabel('max over |dz| <= 0.6 mm (dB)');
legend('1 deg step', '3 deg step', '7.5 deg step', '15 deg step', 'Location', 'northeast');
title('11 angles: transmit replica at \lambda/\Delta(sin\theta)');
dl = Na*vel/PRF/lam;
subplot(1,2,2); plot(dl, res3(:,1), 'k--', dl, res3(:,2), 'o-', dl, res3(:,3), 's-', dl, res3(:,5), 'd-', 'LineWidth', 1.2);
grid on; xlabel('axial displacement v N/PRF over the sequence (\lambda)'); ylabel('coherent loss (dB)');
legend('phase-ramp law', 'sequential: at target', 'sequential: peak', 'alternating: peak', 'Location', 'southwest');
title('motion between angles');
