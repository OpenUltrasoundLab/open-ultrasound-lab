%% m03_lti_eigenfunctions.m
% Purpose : show that sinusoids pass through an LTI filter at the same
%           frequency, scaled by |H| and shifted by angle(H); measure the
%           in-phase and quadrature gains by least squares and compare them
%           with Re H and -Im H; show a complex phasor coming out as H times
%           itself; show that a squarer (nonlinear) creates new frequencies.
% Change  : the filters in hs, the test frequencies ws.
% Expect  : measured in-phase gain = Re H and quadrature gain = -Im H to
%           ~1e-15; for [1/2 1/2] at w = pi/2: 0.5 and 0.5 (H = 0.5 - 0.5j,
%           -45 deg). The squarer output has lines only at 0 and 2*w0.
clear; close all; rng(1);

hs    = {[0.5 0.5], ones(1,3)/3, [1 -1], 0.2*0.8.^(0:79)};
names = {'2-pt average', '3-pt average', 'first difference', 'one-pole a=0.8 (80 taps)'};
ws    = [0.1 0.25 0.5 0.75] * pi;          % rad/sample
n     = 0:399;                              % sample index
Hof   = @(h, w) sum(h .* exp(-1j * w * (0:numel(h)-1)));

%% Cosine in: fit y = a cos(wn) + b sin(wn) on the steady-state part
fprintf('%-26s %6s %9s %9s %9s %9s %8s %8s\n', 'filter', 'w/pi', 'Re H', 'a(fit)', '-Im H', 'b(fit)', '|H|', 'angH');
maxerr = 0;
for i = 1:numel(hs)
    h = hs{i}; M = numel(h);
    for w = ws
        x = cos(w * n);
        y = filter(h, 1, x);                % causal FIR, zero initial state
        keep = n >= M;                      % drop the start-up transient
        A = [cos(w * n(keep)).', sin(w * n(keep)).'];
        ab = A \ y(keep).';
        H = Hof(h, w);
        maxerr = max([maxerr, abs(ab(1) - real(H)), abs(ab(2) + imag(H))]);
        fprintf('%-26s %6.2f %9.4f %9.4f %9.4f %9.4f %8.4f %7.1f\n', names{i}, w/pi, ...
                real(H), ab(1), -imag(H), ab(2), abs(H), angle(H)*180/pi);
    end
end
fprintf('max |fit - prediction| over all cases: %.1e\n', maxerr);

%% Complex phasor in: y./x is the constant H (the eigenvalue)
h = [0.5 0.5]; w = pi/2;
x = exp(1j * w * n); y = filter(h, 1, x);
ratio = y(3:end) ./ x(3:end);
fprintf('phasor in, [1/2 1/2], w = pi/2: mean(y/x) = %s, std(y/x) = %.1e, H = %s\n', ...
        num2str(mean(ratio), 4), std(ratio), num2str(Hof(h, w), 4));
yb = filter(h, 1, cos(w*(-1:6)));          % x(-1) = 0 is included, so y(0..6) is exact
fprintf('  cos(pi n/2) in, n = 0..6: y = %s\n', mat2str(round(yb(2:end)*1e4)/1e4));

%% Nonlinear system: squarer creates DC and 2*w0
N = 200; w0 = 0.3*pi;                       % w0 on DFT bin 30 of 200
xs = cos(w0 * (0:N-1));
S1 = abs(fft(xs)) / N; S2 = abs(fft(xs.^2)) / N;
fprintf('cos(w0 n), w0 = 0.3 pi: |DFT|/N at bins 0, 30, 60: input %.3f %.3f %.3f, squared %.3f %.3f %.3f\n', ...
        S1(1), S1(31), S1(61), S2(1), S2(31), S2(61));
fprintf('  (squared: 1/2 at DC and 1/4 at +-2w0, i.e. cos^2 = 1/2 + 1/2 cos 2w0 n)\n');

%% Figures
wg = linspace(0, pi, 400);
figure('Name', 'Re H and Im H');
subplot(2,1,1);
H1 = arrayfun(@(v) Hof([0.5 0.5], v), wg); H3 = arrayfun(@(v) Hof([1 -1], v), wg);
plot(wg/pi, real(H1), 'b', wg/pi, -imag(H1), 'b--', wg/pi, real(H3), 'r', wg/pi, -imag(H3), 'r--');
xlabel('\omega/\pi'); ylabel('gain');
title('Re H (solid) and -Im H (dashed): 2-pt average (b), first difference (r)'); grid on;
subplot(2,1,2);
h = [0.5 0.5]; w = pi/2; nn = 0:12; tt = linspace(0, 12, 400); H = Hof(h, w);
plot(tt, cos(w*tt), 'Color', [0.5 0.5 0.5]); hold on;
plot(tt, real(H)*cos(w*tt), 'b--', tt, -imag(H)*sin(w*tt), 'r--', tt, abs(H)*cos(w*tt + angle(H)), 'k');
stem(nn, abs(H)*cos(w*nn + angle(H)), 'k', 'filled'); hold off;
xlabel('n (samples)'); title('out (k) = Re H cos (b) - Im H sin (r); input gray');
