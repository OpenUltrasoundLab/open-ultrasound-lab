%% m09_pole_zero_response.m
% Module 09 - The response of a filter read off its pole-zero map.
%  1) |H| = |g| prod|e^jw - z_k| / prod|e^jw - p_k|: geometry vs polynomial evaluation
%  2) A pole at radius r: -3 dB width ~ 2(1-r) rad/sample, peak group delay r/(1-r),
%     ringing time -1/ln(r) ~ 1/(1-r) samples
%  3) Linear-phase FIR: zeros on the unit circle or in reciprocal pairs (z, 1/z*)
% What to change: r_list, theta, the FIR length and cutoff.
% Expected output: geometric and polynomial |H| agree to ~1e-15; for r = 0.95 the
% width is ~0.102 rad/sample (2(1-r) = 0.1) and the peak delay ~18.5 samples for the
% pole pair (19 for a single pole); FIR zeros print as reciprocal pairs.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.

clear; close all; clc;
rng(1);

%% Helpers (anonymous only)
% frequency response and group delay of b/a on a grid w (rad/sample), via
% tau = Re{ sum(n c_n e^-jwn) / sum(c_n e^-jwn) } for numerator minus denominator
freqresp = @(c, w) exp(-1j * w(:) * (0:numel(c)-1)) * c(:);
dfreq    = @(c, w) exp(-1j * w(:) * (0:numel(c)-1)) * ((0:numel(c)-1)' .* c(:));
grpdel   = @(b, a, w) real(dfreq(b, w) ./ freqresp(b, w)) - real(dfreq(a, w) ./ freqresp(a, w));

%% 1) Resonator: geometry vs polynomial
r = 0.95; theta = pi/4;                       % pole radius and angle
p = [r*exp(1j*theta); r*exp(-1j*theta)];      % conjugate pole pair, no zeros
a = real(poly(p)); b = 1;
w = linspace(0, pi, 4001)';
Hpoly = abs(freqresp(b, w) ./ freqresp(a, w));
Hgeo  = 1 ./ (abs(exp(1j*w) - p(1)) .* abs(exp(1j*w) - p(2)));
fprintf('Resonator r = %.2f, theta = pi/4\n', r);
fprintf('  max |Hgeo - Hpoly| / max|H|          = %.2e\n', max(abs(Hgeo - Hpoly)) / max(Hpoly));
fprintf('  |H| at DC = %.3f, at theta = %.2f, peak/DC = %.1f dB\n', Hpoly(1), 1/(abs(exp(1j*theta)-p(1))*abs(exp(1j*theta)-p(2))), 20*log10(max(Hpoly)/Hpoly(1)));

%% 2) Width, delay and ringing versus r
r_list = [0.8 0.9 0.95 0.99];
fprintf('\n   r    -3dB width  2(1-r)   peak tau(pair)  r/(1-r)   -1/ln r\n');
for r = r_list
    p = [r*exp(1j*theta); r*exp(-1j*theta)]; a = real(poly(p));
    H = abs(freqresp(1, w) ./ freqresp(a, w)); [Hm, im] = max(H);
    above = find(H >= Hm/sqrt(2)); width = w(max(above)) - w(min(above));
    tau = grpdel(1, a, w);
    fprintf('  %.2f   %8.4f   %7.4f   %10.2f    %8.2f   %7.2f\n', r, width, 2*(1-r), tau(im), r/(1-r), -1/log(r));
end

%% 3) Linear-phase FIR zeros (21 taps, Kaiser beta = 5, cutoff 0.42 pi)
M = 20; wc = 0.42*pi; beta = 5;
k = 0:40; I0 = @(x) sum(((x(:)/2).^(2*k)) ./ (factorial(k).^2), 2);   % modified Bessel I0 by series
n = (0:M)'; m = n - M/2;
hd = wc/pi * ones(M+1, 1); nz = m ~= 0; hd(nz) = sin(wc*m(nz)) ./ (pi*m(nz));
win = I0(beta*sqrt(1 - (2*n/M - 1).^2)) / I0(beta);
h = hd .* win;
z = roots(h);
fprintf('\nLinear-phase FIR: %d taps, %d zeros\n', M+1, numel(z));
fprintf('  zeros on the unit circle (| |z| - 1 | < 1e-6): %d\n', sum(abs(abs(z) - 1) < 1e-6));
zin = z(abs(z) < 1 - 1e-6); zout = z(abs(z) > 1 + 1e-6);
for i = 1:numel(zin)
    [dmin, j] = min(abs(zout - 1/conj(zin(i))));
    fprintf('  |z| = %.4f  paired with |z| = %.4f   product = %.6f (match err %.1e)\n', abs(zin(i)), abs(zout(j)), abs(zin(i))*abs(zout(j)), dmin);
end
fprintf('  symmetry check max|h - flip(h)| = %.1e, delay M/2 = %d samples\n', max(abs(h - flipud(h))), M/2);

%% Figures
figure('Name', 'Pole-zero geometry');
r = 0.95; p = [r*exp(1j*theta); r*exp(-1j*theta)]; a = real(poly(p));
subplot(2,2,1); tt = linspace(0, 2*pi, 400); plot(cos(tt), sin(tt), 'k-'); hold on;
plot(real(p), imag(p), 'rx', 'MarkerSize', 10, 'LineWidth', 2); plot(cos(theta), sin(theta), 'o', 'MarkerFaceColor', [0.8 0.6 0]);
plot([real(p(1)) cos(theta)], [imag(p(1)) sin(theta)], '-', 'Color', [0.8 0.6 0]); plot([real(p(2)) cos(theta)], [imag(p(2)) sin(theta)], '-', 'Color', [0.8 0.6 0]);
axis image; axis([-1.3 1.3 -1.3 1.3]); grid on; title('Resonator poles, r = 0.95'); xlabel('Re z'); ylabel('Im z');
subplot(2,2,2); plot(w/pi, 20*log10(abs(freqresp(1, w) ./ freqresp(a, w))), 'LineWidth', 1.5); grid on;
xlabel('\omega / \pi'); ylabel('|H| (dB)'); title('Magnitude');
subplot(2,2,3); hold on; cols = lines(4);
for i = 1:numel(r_list)
    rr = r_list(i); pp = [rr*exp(1j*theta); rr*exp(-1j*theta)]; plot(w/pi, grpdel(1, real(poly(pp)), w), 'Color', cols(i,:), 'LineWidth', 1.3);
end
grid on; xlabel('\omega / \pi'); ylabel('\tau_g (samples)'); title('Group delay, r = 0.8 ... 0.99'); ylim([-5 100]);
subplot(2,2,4); plot(cos(tt), sin(tt), 'k-'); hold on; plot(real(z), imag(z), 'bo', 'LineWidth', 1.5);
axis image; axis([-2.5 2.5 -2.5 2.5]); grid on; title('Zeros of a 21-tap linear-phase FIR'); xlabel('Re z'); ylabel('Im z');
