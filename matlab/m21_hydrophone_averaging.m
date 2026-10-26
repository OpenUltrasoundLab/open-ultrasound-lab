%% m21_hydrophone_averaging.m
% Hydrophone spatial averaging at a focus and band-limited measurement of a steepened wave.
% What it shows
%   1) a hydrophone of radius b reads the disc average of p(r) = jinc(r/(lambda F#));
%      loss ~ (pi b / lambda F#)^2 / 16 for small b; harmonics (~p^2) suffer more
%   2) a worked case (5 MHz, F# 2, 0.4 mm) and the same hydrophone at 15 MHz
%   3) a plane Fubini waveform (sigma = 0.8) read through a first-order hydrophone roll-off
% What to change: frequency, F#, diameter, sigma, roll-off frequencies.
% Expected output: -0.215 dB at b/lambdaF# = 0.2, -1.35 dB at 0.5; 0.4 mm at 5 MHz F#2: -0.57 dB;
%   at 15 MHz about -5.2 dB; 3f0 roll-off: p+ about 9 % low, p- about 4 % low.
% Units: mm, MHz.  Base MATLAB / Octave (besselj is part of base MATLAB).

clear; close all; rng(1);
c = 1.54;
jinc = @(x) (abs(x) < 1e-9) + (abs(x) >= 1e-9) .* 2 .* besselj(1, pi*x) ./ (pi*x + (abs(x) < 1e-9));
nr = 40; nt = 64;                                % polar quadrature over the disc
rq = sqrt(((1:nr) - 0.5)/nr); tq = ((1:nt) - 0.5)*2*pi/nt;
[RQ, TQ] = meshgrid(rq, tq);
diskAvg = @(prof, x0, bb) mean(mean(prof(sqrt((x0 + bb*RQ.*cos(TQ)).^2 + (bb*RQ.*sin(TQ)).^2))));

%% 1) loss table versus normalized radius
fprintf('b/(lambda F#)   exact (dB)   small-b law (dB)   2nd harmonic (dB)\n');
bs = [0.05 0.1 0.2 0.3 0.4 0.5];
for bb = bs
    e1 = 20*log10(diskAvg(jinc, 0, bb));
    ap = 20*log10(1 - (pi*bb)^2/16);
    e2 = 20*log10(diskAvg(@(r) jinc(r).^2, 0, bb));
    fprintf('   %4.2f          %7.3f        %7.3f            %7.3f\n', bb, e1, ap, e2);
end

%% 2) worked cases
cases = [5 2 0.4; 15 2 0.4; 7.5 1.5 0.2; 10 2 0.2];
fprintf('\nf (MHz)  F#   d (mm)  lambdaF# (mm)  b/lambdaF#  loss (dB)  pressure  intensity\n');
for k = 1:size(cases, 1)
    f = cases(k,1); Fn = cases(k,2); d = cases(k,3); lF = c/f*Fn; bb = d/2/lF;
    m = diskAvg(jinc, 0, bb);
    fprintf('%5.1f   %4.1f  %5.2f     %6.3f        %5.3f     %6.2f    %5.1f %%   %5.1f %%\n', f, Fn, d, lF, bb, 20*log10(m), 100*(1-m), 100*(1-m^2));
end

% measured profile for the first case
f = 5; Fn = 2; d = 0.4; lF = c/f*Fn; bb = d/2/lF;
u = linspace(-2, 2, 201); tr = jinc(abs(u)); ms = zeros(size(u));
for i = 1:numel(u), ms(i) = diskAvg(jinc, u(i), bb); end
ml = 101:(100 + find(tr(101:end) < 0, 1));                                % main lobe (positive half)
hw = @(p) 2*interp1(p(ml) - p(101)/2, u(ml), 0);                         % -6 dB full width
fprintf('\n-6 dB width true %.3f mm, measured %.3f mm (1.41 lambda F# = %.3f mm)\n', hw(tr)*lF, hw(ms)*lF, 1.4102*lF);
figure(1); clf;
subplot(2,1,1); plot(u*lF, tr, 'k', u*lF, ms, 'b', 'LineWidth', 1.5); grid on;
xlabel('x (mm)'); ylabel('p (norm.)'); title('Focal profile, 5 MHz, F# 2: true (black) and 0.4 mm hydrophone (blue)');
bb2 = linspace(0, 0.6, 31); L1 = zeros(size(bb2)); L2 = L1;
for i = 1:numel(bb2), L1(i) = 20*log10(diskAvg(jinc, 0, bb2(i))); L2(i) = 20*log10(diskAvg(@(r) jinc(r).^2, 0, bb2(i))); end
subplot(2,1,2); plot(bb2, L1, 'b', bb2, 20*log10(1 - (pi*bb2).^2/16), 'k--', bb2, L2, 'r', 'LineWidth', 1.5); grid on;
xlabel('b / (\lambda F#)'); ylabel('peak loss (dB)'); legend('fundamental', 'small-b law', '2nd harmonic', 'Location', 'southwest');

%% 3) band-limited measurement of a Fubini waveform
sig = 0.8; NH = 60; tt = (0:4095)/4096;                                  % one period
Bn = zeros(1, NH); for n = 1:NH, Bn(n) = 2*besselj(n, n*sig)/(n*sig); end
pt = Bn * sin(2*pi*(1:NH)'*tt);
fprintf('\nFubini sigma = %.2f: B2/B1 = %.3f (%.1f dB)\n', sig, Bn(2)/Bn(1), 20*log10(Bn(2)/Bn(1)));
figure(2); clf; plot(tt, pt, 'k', 'LineWidth', 1.6); hold on;
for fcut = [3 5 10]
    nn = (1:NH)'; mag = 1./sqrt(1 + (nn/fcut).^2); ph = -atan(nn/fcut);
    pm = (Bn' .* mag)' * sin(2*pi*nn*tt + repmat(ph, 1, numel(tt)));
    fprintf('roll-off %2d f0: p+ %5.1f %% low, p- %5.1f %% low\n', fcut, 100*(1 - max(pm)/max(pt)), 100*(1 - min(pm)/min(pt)));
    plot(tt, pm, 'LineWidth', 1.1);
end
xlabel('t / period'); ylabel('p / p_0'); grid on; title('Fubini waveform, \sigma = 0.8: true and band-limited hydrophone');
legend('true', '3 f_0', '5 f_0', '10 f_0', 'Location', 'southwest');
