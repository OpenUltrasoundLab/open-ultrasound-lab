%% m22_hydrophone_spatial_averaging.m  -  How a finite hydrophone under-reads a focus, and a correction
% Purpose : average a focal pressure field over a disc of diameter d (the hydrophone face) and
%           report the pressure (MI) and intensity (I_SPTA, I_SPPA) errors versus d/W6, where W6 is
%           the -6 dB beam width; compare with the closed form for a round Gaussian focus; apply
%           the Gaussian-kernel correction that uses only the measured width; repeat for a jinc
%           (uniform circular aperture) focus with the same W6.
% Change  : W6 (mm), aspect (elevation/azimuth width), the list of d/W6 ratios.
% Expect  : round Gaussian, d = W6/4: -0.19 dB in p, -0.38 dB (4.2 %) in I;
%           d = W6/2: -0.74 dB, -1.48 dB (15.7 %); corrected within 0.02 dB at d = W6/2.
rng(1);
W6 = 1.0;                          % mm, azimuth -6 dB width (pressure or intensity: same point)
aspect = 3;                        % elevation width / azimuth width
w = W6/1.6651;                     % mm, 1/e amplitude half-width: p = exp(-r^2/w^2)
ratios = [0.1 0.25 0.5 0.67 0.75 1.0 1.2];
nr = 40; nt = 64;                  % polar quadrature over the disc
discavg = @(fun, b) sum(sum(fun(((1:nr)'-0.5)/nr*b*cos(2*pi*(0:nt-1)/nt), ((1:nr)'-0.5)/nr*b*sin(2*pi*(0:nt-1)/nt)) .* repmat(((1:nr)'-0.5), 1, nt))) / (nt*sum((1:nr)-0.5));
g_round = @(x, y) exp(-(x.^2 + y.^2)/w^2);
g_elong = @(x, y) exp(-x.^2/w^2 - y.^2/(aspect*w)^2);
% jinc focus with the same -6 dB width: p = 2 J1(u)/u, u = k r; -6 dB at u = 2.2151
u6 = fzero(@(u) 2*besselj(1,u)/u - 0.5, 2.2);
kj = u6/(W6/2);
jincf = @(x, y) 2*besselj(1, kj*sqrt(x.^2 + y.^2 + 1e-12))./(kj*sqrt(x.^2 + y.^2 + 1e-12));
fprintf('-6 dB point of 2J1(u)/u at u = %.4f\n', u6);
fprintf(' d/W6 | round: num  closed  (dB p) | I error %% | corrected (dB) | elongated %g:1 (dB p) | jinc (dB p)\n', aspect);
res = zeros(numel(ratios), 5);
for k = 1:numel(ratios)
  b = ratios(k)*W6/2;
  pr = discavg(g_round, b);
  pcf = (w^2/b^2)*(1 - exp(-b^2/w^2));
  % Gaussian-kernel correction: measured width^2 = w^2 + b^2/2 per axis, peak scales by w/w_meas per axis
  cr = pr*(w^2 + b^2/2)/w^2;
  pe = discavg(g_elong, b);
  pj = discavg(jincf, b);
  res(k,:) = [20*log10(pr), 20*log10(pcf), 20*log10(cr), 20*log10(pe), 20*log10(pj)];
  fprintf(' %4.2f |      %6.3f  %6.3f      |  %5.1f    |    %6.3f      |        %6.3f         |  %6.3f\n', ...
    ratios(k), res(k,1), res(k,2), 100*(1 - pr^2), res(k,3), res(k,4), res(k,5));
end
% Example: 0.4 mm hydrophone on a 0.6 mm focus (7.5 MHz, F/2)
b = 0.2; wx = 0.6/1.6651;
p1 = discavg(@(x, y) exp(-(x.^2 + y.^2)/wx^2), b);
p3 = discavg(@(x, y) exp(-x.^2/wx^2 - y.^2/(3*wx)^2), b);
fprintf('\n0.4 mm hydrophone, W6 = 0.6 mm: round %.2f dB (I %.1f %% low), 3:1 elongated %.2f dB (I %.1f %% low)\n', ...
  20*log10(p1), 100*(1-p1^2), 20*log10(p3), 100*(1-p3^2));
%% Measured lateral profile for one case (d = 0.5 W6, round)
b = 0.25*W6; xs = linspace(-1.6*W6, 1.6*W6, 161); meas = zeros(size(xs));
for i = 1:numel(xs)
  meas(i) = discavg(@(x, y) g_round(x + xs(i), y), b);
end
c6 = meas(81)/2; i6 = find(meas(81:end) < c6, 1) + 79;
W6m = 2*interp1(meas(i6:i6+1), xs(i6:i6+1), c6);
fprintf('d = W6/2: measured W6 = %.3f mm (true %.3f), predicted sqrt(W6^2 + 1.6651^2 b^2/2) = %.3f mm\n', W6m, W6, 1.6651*sqrt(w^2 + b^2/2));
%% Figure
figure('Position', [100 100 900 340]);
subplot(1,2,1);
plot(ratios, res(:,1), 'b-o', ratios, 2*res(:,1), 'r-o', ratios, res(:,4), 'b--', ratios, res(:,5), 'k-.', ratios, res(:,3), 'g-s', 'LineWidth', 1.3);
grid on; xlabel('hydrophone diameter / W_{-6}'); ylabel('reading error (dB)');
legend('pressure, round', 'intensity, round', 'pressure, elongated', 'pressure, jinc focus', 'corrected pressure', 'Location', 'southwest');
title('spatial-averaging error');
subplot(1,2,2);
plot(xs, g_round(xs, 0), 'k', xs, meas, 'b', 'LineWidth', 1.5); grid on;
xlabel('x (mm)'); ylabel('p (norm.)'); legend('true', 'measured, d = W_{-6}/2');
title(sprintf('W_{-6} = %.2f mm reads as %.3f mm', W6, W6m));
