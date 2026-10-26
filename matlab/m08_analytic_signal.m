%% m08_analytic_signal.m  -  FFT-based analytic signal of an RF echo
% Purpose : build the analytic signal x_a = x + j*H{x} with an FFT mask
%           (no hilbert()), then take envelope |x_a|, instantaneous phase
%           and instantaneous frequency of a Gaussian-modulated RF echo.
%           Compares the envelope with the true Gaussian envelope, checks the
%           mean-frequency identity, and shows three failure modes:
%           Bedrosian (low f0, very broadband), circular edge wrap, and
%           instantaneous-frequency spikes where the envelope is small.
% Change  : f0, bw (fractional -6 dB bandwidth), t0 (echo position), noise.
% Expect  : envelope rms error ~0.00 % at 5 MHz / 60 %; ~2.5 % of the lobe
%           below 0 Hz at 2 MHz / 120 % (envelope error of a few %);
%           A^2-weighted mean f_i = spectral centroid; truncated echo leaks
%           energy to the far end of the record.
clear; close all; rng(1);

%% Parameters
fs  = 40;            % sampling rate (MHz = samples per us)
N   = 1024;          % record length (samples) -> 25.6 us
t   = (0:N-1)'/fs;   % time (us)
f0  = 5;             % center frequency (MHz)
bw  = 0.60;          % -6 dB fractional bandwidth of the amplitude spectrum
t0  = 12.8;          % echo center (us)

sigf = bw*f0/(2*sqrt(2*log(2)));   % std of Gaussian amplitude spectrum (MHz)
sigt = 1/(2*pi*sigf);              % std of Gaussian envelope (us)

% FFT analytic-signal mask: DC and Nyquist x1, positive bins x2, negative x0
h = zeros(N,1); h(1) = 1; h(2:N/2) = 2; h(N/2+1) = 1;
analytic = @(x) ifft(fft(x).*h);   % works column-wise on an N-row matrix

% Gaussian-modulated echo and its true envelope
echo = @(f0_, sigt_, tc) exp(-(t-tc).^2/(2*sigt_^2)).*cos(2*pi*f0_*(t-tc));
genv = @(sigt_, tc) exp(-(t-tc).^2/(2*sigt_^2));

%% 1. Envelope, phase, instantaneous frequency of a clean echo
x   = echo(f0, sigt, t0);
xa  = analytic(x);
A   = abs(xa);                          % envelope
phi = unwrap(angle(xa));                % instantaneous phase (rad)
fi  = [diff(phi); NaN]*fs/(2*pi);       % instantaneous frequency (MHz)
tr  = genv(sigt, t0);

m   = tr > 0.01;
errEnv = 100*sqrt(sum((A(m)-tr(m)).^2)/sum(tr(m).^2));
fprintf('Re{x_a} - x max abs        = %.2e (should be ~1e-16)\n', max(abs(real(xa)-x)));
fprintf('Envelope rms error (5 MHz, 60 %%)  = %.4f %%\n', errEnv);

% Mean-frequency identity: A^2-weighted f_i  vs  centroid of |X_a|^2
w   = A(1:end-1).*A(2:end);
fiW = sum(w.*fi(1:end-1))/sum(w);
Xa  = fft(xa); fk = (0:N-1)'*fs/N; pos = 1:N/2;
fC  = sum(fk(pos).*abs(Xa(pos)).^2)/sum(abs(Xa(pos)).^2);
fprintf('A^2-weighted mean f_i       = %.4f MHz\n', fiW);
fprintf('Centroid of |X_a(f)|^2      = %.4f MHz\n', fC);
fprintf('Theory (Gaussian, no clip)  = %.4f MHz\n', f0);

%% 2. Bedrosian: low f0, very broadband pulse
f0b = 2; bwb = 1.2; sigfb = bwb*f0b/(2*sqrt(2*log(2))); sigtb = 1/(2*pi*sigfb);
xb  = echo(f0b, sigtb, t0); Ab = abs(analytic(xb)); trb = genv(sigtb, t0);
mb  = trb > 0.01;
errB = 100*sqrt(sum((Ab(mb)-trb(mb)).^2)/sum(trb(mb).^2));
spill = 0.5*erfc(f0b/sigfb/sqrt(2));      % Phi(-f0/sigma_f)
fprintf('\n2 MHz, 120 %%: sigma_f = %.3f MHz, f0/sigma_f = %.2f\n', sigfb, f0b/sigfb);
fprintf('  lobe fraction below 0 Hz   = %.2f %%\n', 100*spill);
fprintf('  envelope rms error         = %.2f %%\n', errB);

%% 3. Circular edge artifact: echo cut by the start of the record
tE  = 0.15;                              % echo center near the record start (us)
xe  = echo(f0, sigt, tE); Ae = abs(analytic(xe));
xe0 = [xe; zeros(N,1)];                  % same echo, zero-padded to 2N
h2  = zeros(2*N,1); h2(1) = 1; h2(2:N) = 2; h2(N+1) = 1;
Ae0 = abs(ifft(fft(xe0).*h2)); Ae0 = Ae0(1:N);
tail = t > t(end) - 0.5;                 % last 0.5 us of the record
fprintf('\nTruncated echo at %.2f us: max |x_a| in last 0.5 us = %.3f (no pad)\n', tE, max(Ae(tail)));
fprintf('                              with 2x zero pad   = %.3f\n', max(Ae0(tail)));

%% 4. Noise: instantaneous frequency is erratic where the envelope is small
xn  = x + 0.05*randn(N,1);
xan = analytic(xn); An = abs(xan);
fin = angle(xan(2:end).*conj(xan(1:end-1)))*fs/(2*pi);   % no unwrap needed
wn  = An(1:end-1).*An(2:end);
win = abs(t(1:end-1) - t0) < 3*sigt;                     % echo window, +-3 sigma_t
strong = An(1:end-1) > 0.5*max(An);
weak   = An(1:end-1) < 0.1*max(An);
fprintf('\nNoise 5 %% of peak:\n');
fprintf('  std f_i where A > 50 %% max = %.2f MHz\n', std(fin(strong)));
fprintf('  std f_i where A < 10 %% max = %.2f MHz\n', std(fin(weak)));
fprintf('  A^2-weighted mean f_i, echo window = %.3f MHz (clean: %.3f)\n', sum(wn(win).*fin(win))/sum(wn(win)), fiW);
fprintf('  unweighted mean f_i, echo window   = %.3f MHz\n', mean(fin(win)));

%% Figures
figure('Position', [100 100 900 650]);
zw = abs(t - t0) < 6*sigt;
subplot(2,2,1);
plot(t(zw), x(zw), 'b', t(zw), imag(xa(zw)), 'r', t(zw), A(zw), 'k', t(zw), -A(zw), 'k', ...
     t(zw), tr(zw), 'g--');
xlabel('time (\mus)'); ylabel('amplitude');
title(sprintf('x (b), Hilbert (r), |x_a| (k), Gaussian (g--); err %.3f %%', errEnv));
subplot(2,2,2);
fz = zw & A > 0.1*max(A);
plot(t(fz), fi(fz), 'b.-'); hold on; plot(t(fz), f0+0*t(fz), 'k--');
xlabel('time (\mus)'); ylabel('f_i (MHz)'); title('instantaneous frequency (A > 10 % max)');
ylim([f0-1 f0+1]);
subplot(2,2,3);
plot(t, 20*log10(Ae+1e-12), 'r', t, 20*log10(Ae0+1e-12), 'k--'); ylim([-80 5]);
xlabel('time (\mus)'); ylabel('|x_a| (dB)');
title('truncated echo: no pad (r) vs 2x zero pad (k--)'); legend('no pad', 'zero pad');
subplot(2,2,4);
plot(t(zw), fin(min(find(zw), N-1)), 'b.'); hold on;
plot(t(zw), An(zw)/max(An)*4 + f0 - 2, 'k');
xlabel('time (\mus)'); ylabel('f_i (MHz)'); ylim([-5 15]);
title('noisy echo: f_i (dots) and scaled envelope (k)');
