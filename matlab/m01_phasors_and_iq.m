%% m01_phasors_and_iq.m  --  rotating phasors, the imaginary part, and I/Q of an RF echo
% Purpose : sections 01.4-01.6 by computer.
%   (1) e^{j 2 pi f t}: angle at a given time, the sign of f flips sin but not cos
%   (2) cos = (e^{+j} + e^{-j})/2: the imaginary parts cancel to machine precision
%   (3) the imaginary part (Q) of a narrowband signal = the real signal a quarter period earlier
%   (4) I/Q by hand at fs = 4 f0 (multiply by 1, -j, -1, j, average, x2)
%   (5) I/Q demodulation of a Gaussian RF echo: mix with e^{-j 2 pi fd t}, low-pass, x2;
%       phase at the peak = -2 pi f0 tau, envelope = |IQ|; residual rotation if fd ~= f0
% Change  : tau (echo delay), fd (demodulation frequency), bw (fractional -6 dB bandwidth).
% Expect  : angle at 50 ns = 90 deg; hand I/Q = 0.8660 + 0.5000j = 1 at 30 deg;
%           tau = 13 ns -> phase -23.4 deg; fd = 4.5 MHz -> IQ spins at +0.5 MHz.
% Runtime : < 2 s. Base MATLAB / GNU Octave, no toolboxes.
clear; close all; rng(1);

%% 1. A rotating arrow
f0 = 5e6;                         % carrier (Hz) = turns per second
t1 = 50e-9;                       % 50 ns = a quarter of the 200 ns period
p_pos = exp(1j*2*pi*f0*t1);
p_neg = exp(-1j*2*pi*f0*t1);
fprintf('f0 = 5 MHz, period %.0f ns. At t = 50 ns the angle is %.1f deg\n', 1e9/f0, 2*pi*f0*t1*180/pi);
fprintf('  +f: %+.3f %+.3fj ;  -f: %+.3f %+.3fj  (same cos, opposite sin)\n\n', real(p_pos), imag(p_pos), real(p_neg), imag(p_neg));

%% 2. Two counter-rotating halves make a real cosine
t = (0:999)*1e-10;                % 0 .. 100 ns, 0.1 ns steps
x = 0.5*exp(1j*2*pi*f0*t) + 0.5*exp(-1j*2*pi*f0*t);
fprintf('max |imag( e^{+j}/2 + e^{-j}/2 )| = %.1e ; max |real - cos| = %.1e\n', max(abs(imag(x))), max(abs(real(x) - cos(2*pi*f0*t))));
x_unbal = 0.5*exp(1j*2*pi*f0*t) + 0.25*exp(-1j*2*pi*f0*t);   % unbalanced: no longer real
fprintf('unbalanced (0.5, 0.25): max |imag| = %.3f = 0.5 - 0.25\n\n', max(abs(imag(x_unbal))));

%% 3. The imaginary part is the quarter-period-earlier copy
A = 1; phi = 30*pi/180; T0 = 1/f0;
xr = @(tt) A*cos(2*pi*f0*tt + phi);       % the real RF
I = A*cos(phi); Q = A*sin(phi);           % phasor A e^{j phi} = I + jQ
fprintf('A cos(wt + 30 deg): I = A cos(phi) = %.4f, Q = A sin(phi) = %.4f\n', I, Q);
fprintf('  x(0) = %.4f (= I);  x(-T/4) = %.4f (= Q)\n', xr(0), xr(-T0/4));
fprintf('  and x(t) = I cos(wt) + Q cos(wt + 90 deg): error %.1e\n\n', max(abs(xr(t) - (I*cos(2*pi*f0*t) + Q*cos(2*pi*f0*t + pi/2)))));

%% 4. I/Q by hand at fs = 4 f0
n = 0:3;
s  = cos(n*pi/2 + phi);                   % four samples, one per quarter period
lo = exp(-1j*n*pi/2);                     % e^{-j 2 pi f0 t_n} = 1, -j, -1, j
prod_ = s .* lo;
iq = 2*mean(prod_);
fprintf('fs = 4 f0 samples   : %s\n', sprintf('%+.4f  ', s));
fprintf('times 1, -j, -1, j  : %s\n', sprintf('%+.4f%+.4fj  ', [real(prod_); imag(prod_)]));
fprintf('2 x mean = %.4f %+.4fj  -> |IQ| = %.4f, angle %.2f deg\n\n', real(iq), imag(iq), abs(iq), angle(iq)*180/pi);

%% 5. I/Q demodulation of an RF echo
fs  = 40e6;                       % sampling rate (Hz)
N   = 256;                        % 6.4 us record
tt  = (0:N-1)/fs;
t0  = 3.2e-6;                     % reference echo time: 16 whole periods of f0
bw  = 0.6;                        % -6 dB fractional bandwidth (60 %)
sig = sqrt(log(2)/2)/(pi*bw*f0/2);% Gaussian envelope std (s): 125 ns for 60 %
tau = 13e-9;                      % extra delay (s)
rf  = exp(-(tt - t0 - tau).^2/(2*sig^2)) .* cos(2*pi*f0*(tt - t0 - tau));
fk  = [0:N/2-1, -N/2:-1]*fs/N;    % FFT frequency axis (Hz)
lpf = double(abs(fk) < 5e6);      % ideal low-pass, cutoff 5 MHz (halfway to the 2f0 image)
demod = @(fd) 2*ifft(fft(rf .* exp(-1j*2*pi*fd*tt)) .* lpf);
iq0 = demod(f0);
[~, ipk] = max(abs(iq0));
fprintf('Gaussian echo: f0 = 5 MHz, 60 %% bandwidth (sigma = %.0f ns), fs = 40 MHz, tau = %.0f ns\n', sig*1e9, tau*1e9);
fprintf('  fd = f0  : phase at the peak = %.2f deg ; predicted -360 f0 tau = %.2f deg\n', angle(iq0(ipk))*180/pi, -360*f0*tau);
fprintf('             peak |IQ| = %.4f (envelope peak 1)\n', abs(iq0(ipk)));
fd2 = 4.5e6; iq2 = demod(fd2);
win = abs(tt - t0) < 1.5*sig;                         % inside the echo
ph2 = unwrap(angle(iq2(win)));
rate = (ph2(end) - ph2(1)) / (2*pi*(sum(win) - 1)/fs); % turns per second
fprintf('  fd = 4.5 MHz: IQ keeps spinning at %.3f MHz (= f0 - fd); max envelope change %.1e\n', rate/1e6, max(abs(abs(iq2) - abs(iq0))));
% quarter-period trick on the broadband echo (no Hilbert filter), error relative to the analytic signal
for bwq = [0.1 0.6]               % fractional -6 dB bandwidth
    sq  = sqrt(log(2)/2)/(pi*bwq*f0/2);                  % envelope std (s)
    rq  = exp(-(tt - t0).^2/(2*sq^2)) .* cos(2*pi*f0*(tt - t0));
    xaq = ifft(fft(rq) .* (1 + sign(fk)));               % analytic signal (keep +f, x2)
    xqq = interp1(tt, rq, tt - T0/4, 'spline', 0);       % real signal a quarter period earlier
    fprintf('  quarter-period copy vs Hilbert at %2.0f %% bandwidth: rms error %.1f %% of the analytic-signal rms\n', 100*bwq, 100*sqrt(sum(abs(rq + 1j*xqq - xaq).^2)/sum(abs(xaq).^2)));
end

%% 6. Figures
figure('Position', [80 80 1100 380]);
m = tt > 2.6e-6 & tt < 3.8e-6;
subplot(1,3,1);
plot(tt(m)*1e6, rf(m), 'o-', 'Color', [0.12 0.37 0.82], 'MarkerSize', 3); hold on; grid on;
plot(tt(m)*1e6, abs(iq0(m)), 'k--', tt(m)*1e6, -abs(iq0(m)), 'k--');
xlabel('time (\mus)'); ylabel('amplitude'); title('RF echo at 40 MHz and envelope |IQ|');
subplot(1,3,2);
plot(tt(m)*1e6, real(iq0(m)), '-', 'Color', [0.12 0.37 0.82], 'LineWidth', 1.5); hold on; grid on;
plot(tt(m)*1e6, imag(iq0(m)), '-', 'Color', [0.76 0.21 0.44], 'LineWidth', 1.5);
plot(tt(m)*1e6, real(iq2(m)), ':', 'Color', [0.12 0.37 0.82]);
legend('I (fd = f0)', 'Q (fd = f0)', 'I (fd = 4.5 MHz)', 'Location', 'southwest');
xlabel('time (\mus)'); ylabel('amplitude'); title('I and Q: the slow arrow');
subplot(1,3,3);
plot(real(iq0(m)), imag(iq0(m)), '-', 'Color', [0.12 0.37 0.82], 'LineWidth', 3); hold on; grid on;
plot(real(iq2(m)), imag(iq2(m)), '-', 'Color', [0.66 0.43 0.02]);
plot(real(iq0(ipk)), imag(iq0(ipk)), 'ko', 'MarkerSize', 7);
axis equal; axis([-1.1 1.1 -1.1 1.1]); xlabel('I'); ylabel('Q');
legend('fd = f0: one fixed direction', 'fd = 4.5 MHz: loop (residual turn)', 'Location', 'southoutside');
title('IQ plane: phase -23.4 deg at the peak');
