%% m01_phase_and_derivative.m  --  adding phasors, beats, phase unwrapping, d/dt = j omega
% Purpose : sections 01.7-01.9 by computer.
%   (1) equal-frequency sinusoids add as arrows (phasor sum vs brute-force time sum)
%   (2) two scatterers dz apart in depth: |1 + e^{-j 4 pi dz / lambda}| = 2|cos(2 pi dz/lambda)|
%   (3) beats: 5.0 + 5.5 MHz, envelope period 1/df = 2 us
%   (4) unwrap works only while the true phase step per sample is below pi
%   (5) differentiation multiplies each phasor by j omega: +90 deg and +6.02 dB per octave
% Change  : the phasors in (1), df in (3), the acceleration a in (4).
% Expect  : cos + cos(+90) = 1.4142 cos(wt + 45 deg); 1 + 0.5 at 120 deg = 0.8660 at 30 deg;
%           cancellation at dz = lambda/4 = 77 um (5 MHz); beat period 2.00 us;
%           unwrap first fails at sample 32 for a = 0.1 rad/sample^2; derivative gain 6.02 dB per octave;
%           Gaussian pulse spectral peak 5.00 -> 5.46 MHz after d/dt.
% Runtime : < 2 s. Base MATLAB / GNU Octave, no toolboxes.
clear; close all; rng(1);

%% 1. Same frequency: add the arrows
f0 = 5e6; c = 1540; lambda = c/f0;            % 0.308 mm
t  = (0:3999)*1e-10;                           % 0 .. 400 ns
pairs = {[1, 0; 1, 90], [1, 0; 0.5, 120], [1, 0; 1, 120; 1, 240], [1, 0; 1, 180]};
for k = 1:numel(pairs)
    P = pairs{k};                              % rows: [amplitude, phase in deg]
    S = sum(P(:,1) .* exp(1j*P(:,2)*pi/180));  % phasor sum
    xsum = zeros(size(t));
    for r = 1:size(P,1), xsum = xsum + P(r,1)*cos(2*pi*f0*t + P(r,2)*pi/180); end
    err = max(abs(xsum - abs(S)*cos(2*pi*f0*t + angle(S))));
    if abs(S) < 1e-12, S = 0; end             % exact cancellation: the angle is undefined
    fprintf('phasors %-28s -> sum %.4f at %7.2f deg   (time-domain check %.1e)\n', mat2str(P), abs(S), angle(S)*180/pi, err);
end

%% 2. Two scatterers in one resolution cell
dz = [0 0.25 0.5 1] * lambda/2;               % axial separations
amp = abs(1 + exp(-1j*4*pi*dz/lambda));       % round trip: phase 4 pi dz / lambda
fprintf('\nlambda = %.3f mm. Two equal echoes, separation dz:\n', lambda*1e3);
for i = 1:numel(dz)
    fprintf('  dz = %6.1f um  (%4.2f lambda)  ->  |sum| = %.3f\n', dz(i)*1e6, dz(i)/lambda, amp(i));
end

%% 3. Beats
fs = 200e6; tb = (0:1999)/fs;                 % 10 us: whole periods of both tones
df = 0.5e6;
xb = cos(2*pi*f0*tb) + cos(2*pi*(f0 + df)*tb);
fk = [0:999, -1000:-1]*fs/2000;
xa = ifft(fft(xb) .* (1 + sign(fk)));         % analytic signal
env = abs(xa);
mins = find(env(2:end-1) < env(1:end-2) & env(2:end-1) < env(3:end) & env(2:end-1) < 0.2) + 1;
fprintf('\nbeats 5.0 + 5.5 MHz: envelope minima spacing %.3f us (1/df = %.3f us); max %.3f, min %.3f\n', ...
        mean(diff(tb(mins)))*1e6, 1e6/df, max(env(200:end-200)), min(env(200:end-200)));

%% 4. Phase wrapping and unwrapping
a  = 0.1;                                     % rad per sample^2 : phase step grows by a each sample
n  = 0:59;
ph_true = 0.5*a*n.^2;                         % step between samples n-1 and n = a(n - 0.5)
ph_wrap = angle(exp(1j*ph_true));             % what angle() returns, in (-pi, pi]
ph_unw  = unwrap(ph_wrap);
bad = find(abs(ph_unw - ph_true) > 1e-9, 1);
fprintf('\nunwrap: a = %.2f rad/sample^2, first wrong sample n = %d (true step there %.3f rad > pi)\n', a, n(bad), ph_true(bad) - ph_true(bad-1));
fprintf('        at 5 MHz a step of pi is a round-trip displacement of lambda/4 = %.1f um per sample\n', lambda/4*1e6);
fprintf('        wrapped values of 0, 2, 4, 6 rad: %s\n', sprintf('%.4f  ', angle(exp(1j*[0 2 4 6]))));

%% 5. d/dt is multiplication by j omega
fs2 = 400e6; N = 4096; td = ((0:N-1) - N/2)/fs2;
sig = 100e-9;
g  = exp(-td.^2/(2*sig^2)) .* cos(2*pi*f0*td);           % Gaussian RF pulse
dg_exact = (-td/sig^2).*exp(-td.^2/(2*sig^2)).*cos(2*pi*f0*td) - 2*pi*f0*exp(-td.^2/(2*sig^2)).*sin(2*pi*f0*td);
fk2 = [0:N/2-1, -N/2:-1]*fs2/N;
dg_fft = real(ifft(1j*2*pi*fk2 .* fft(g)));               % multiply each phasor by j 2 pi f
fprintf('\nFFT derivative vs exact: max error %.1e of max %.3e /s\n', max(abs(dg_fft - dg_exact)), max(abs(dg_exact)));
fprintf('|H| = 2 pi f: 5 MHz -> %.3e /s, 10 MHz -> %.3e /s, ratio %.2f dB (one octave)\n', 2*pi*5e6, 2*pi*10e6, 20*log10(2));
tc = (0:799)/fs2;                                        % exactly 10 periods of f0
fc = [0:399, -400:-1]*fs2/800;
d1 = real(ifft(1j*2*pi*fc .* fft(cos(2*pi*f0*tc))));     % derivative of a pure cosine
fprintf('d/dt cos: phase of H at f0 = %.1f deg ; |d/dt cos(wt) - w cos(wt + 90 deg)| / w = %.1e\n', ...
        angle(1j*2*pi*f0)*180/pi, max(abs(d1 - 2*pi*f0*cos(2*pi*f0*tc + pi/2)))/(2*pi*f0));
G = abs(fft(g)); D = abs(fft(dg_exact));
Np = 2^17; fkp = (0:Np/2-1)*fs2/Np;                     % zero-padded grid (3 kHz bins): peaks not bin-quantized
Gp = abs(fft(g, Np)); Dp = abs(fft(dg_exact, Np));
[~, ig] = max(Gp(1:Np/2)); [~, id] = max(Dp(1:Np/2));
sigf = 1/(2*pi*sig);                                     % spectral std of the Gaussian envelope, Hz
fprintf('spectral peak: pulse %.2f MHz, its derivative %.2f MHz (the jw weighting raises it)\n', fkp(ig)/1e6, fkp(id)/1e6);
fprintf('               predicted (f0 + sqrt(f0^2 + 4 sigf^2))/2 = %.2f MHz, sigf = %.2f MHz\n', (f0 + sqrt(f0^2 + 4*sigf^2))/2/1e6, sigf/1e6);

%% 6. Figures
figure('Position', [80 80 1100 380]);
subplot(1,3,1);
plot(tb*1e6, xb, 'Color', [0.12 0.37 0.82]); hold on; grid on;
plot(tb*1e6, env, 'k', 'LineWidth', 1.5);
xlim([0 6]); xlabel('time (\mus)'); ylabel('amplitude');
title('beats: 5.0 + 5.5 MHz, envelope period 2 \mus');
subplot(1,3,2);
plot(n, ph_true, 'k-', 'LineWidth', 1.5); hold on; grid on;
plot(n, ph_wrap, 'o', 'Color', [0.76 0.21 0.44], 'MarkerSize', 4);
plot(n, ph_unw, '-', 'Color', [0.12 0.37 0.82], 'LineWidth', 1.5);
plot([n(bad) n(bad)], [-5 180], ':', 'Color', [0.66 0.43 0.02]);
xlabel('sample n'); ylabel('phase (rad)'); legend('true', 'wrapped', 'unwrapped', 'step > \pi', 'Location', 'northwest');
title('unwrap fails once the step exceeds \pi');
subplot(1,3,3);
ff = fk2(1:N/2)/1e6;
plot(ff, G(1:N/2)/max(G), 'Color', [0.12 0.37 0.82], 'LineWidth', 1.5); hold on; grid on;
plot(ff, D(1:N/2)/max(D), 'Color', [0.76 0.21 0.44], 'LineWidth', 1.5);
xlim([0 12]); xlabel('frequency (MHz)'); ylabel('normalized |spectrum|');
legend('pulse', 'd/dt pulse', 'Location', 'northeast'); title('d/dt weights each line by 2\pi f');
