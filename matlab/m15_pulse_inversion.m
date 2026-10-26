%% m15_pulse_inversion.m
% Purpose : pulse inversion (PI) and amplitude modulation (AM) as Volterra-order filters.
%           Tissue: the transmitted pulse distorts by lossless plane-wave propagation (simple-wave
%           solution, each sample of the waveform travels at c0 + beta u), then scatters linearly.
%           "Bubble-like": a memoryless polynomial echo y = x + a2 x^2 + a3 x^3 (illustrative only).
%           PI  : e(+s) + e(-s)          -> odd orders cancel, even orders double
%           AM  : e(s) - 2 e(s/2)        -> linear cancels, quadratic x 1/2, cubic x 3/4
%           Motion between pulses leaves 2|sin(pi f0 delta)| of the fundamental, delta = 2 v / (c PRF).
% Change  : f0, FBW, sigma (distortion = x / xbar), a2, a3, velocity and PRF.
% Expect  : tissue, sigma = 0.3: single echo 2f0/f0 band ratio about -18 dB; PI removes the fundamental
%           band to the level of the even-order spill (< -30 dB re single), 2f0 band +5.8 dB (ideal +6.0);
%           motion residual -26.2 dB at 10 mm/s, PRF 5 kHz, 3 MHz (formula and simulation agree).
%           Band separation: Gaussian f0 and 2f0 spectra meet at -20 dB when FBW = 45 %.
% Runtime : under 1 s.
clear; close all; rng(1);

%% Pulse and helpers
c = 1540; f0 = 3e6; FBW = 0.5; fs = 200e6; N = 4096;
t = ((0:N-1) - N/2)/fs; sf = FBW*f0/(2*sqrt(2*log(2))); st = 1/(2*pi*sf);
s = exp(-t.^2/(2*st^2)).*sin(2*pi*f0*t);                % unit-amplitude transmit pulse
fk = (0:N-1)*fs/N;
bandE = @(x, f1, f2) sum(abs(fft(x)).^2 .* (fk >= f1 & fk < f2));
Ef = @(x) bandE(x, 0.5*f0, 1.5*f0); Eh = @(x) bandE(x, 1.5*f0, 2.5*f0);
dB10 = @(r) 10*log10(r);
% lossless simple-wave distortion: a sample of value A*s(t) arrives at t - sigma*A*s(t)/(2 pi f0)
% (sigma is defined for amplitude 1; amplitude A scales the distortion by A). Valid before the shock.
prop = @(A, sig) interp1(t - sig*A*s/(2*pi*f0), A*s, t, 'linear', 0);

%% Tissue: PI and AM
sig = 0.3;
eP = prop(1, sig); eN = prop(-1, sig); eH = prop(0.5, sig);
pi_out = eP + eN; am_out = eP - 2*eH;
fprintf('Tissue, f0 = %.1f MHz, FBW = %.0f %%, sigma = %.2f\n', f0/1e6, FBW*100, sig);
fprintf(' single echo: 2f0 band re f0 band          %6.1f dB\n', dB10(Eh(eP)/Ef(eP)));
fprintf(' PI : f0 band re single f0 band             %6.1f dB\n', dB10(Ef(pi_out)/Ef(eP)));
fprintf(' PI : 2f0 band re single 2f0 band           %6.1f dB  (ideal +6.0)\n', dB10(Eh(pi_out)/Eh(eP)));
fprintf(' AM : f0 band re single f0 band             %6.1f dB\n', dB10(Ef(am_out)/Ef(eP)));
fprintf(' AM : 2f0 band re single 2f0 band           %6.1f dB  (quadratic retained x 1/2: -6.0)\n', dB10(Eh(am_out)/Eh(eP)));

%% Bubble-like polynomial echo
a2 = 0.3; a3 = 0.2; ech = @(x) x + a2*x.^2 + a3*x.^3;
bP = ech(s); bN = ech(-s); bH = ech(0.5*s);
fprintf('\nPolynomial echo (a2 = %.1f, a3 = %.1f): PI keeps even orders, AM also keeps cubic terms\n', a2, a3);
fprintf(' PI : f0 band re single f0 band             %6.1f dB\n', dB10(Ef(bP + bN)/Ef(bP)));
fprintf(' AM : f0 band re single f0 band             %6.1f dB  (nonlinear fundamental from x^3)\n', dB10(Ef(bP - 2*bH)/Ef(bP)));
fprintf(' AM : expected cubic-fundamental level      %6.1f dB\n', dB10(Ef(0.75*a3*s.^3)/Ef(bP)));

%% Motion between pulses
fprintf('\nMotion residual of the fundamental after PI (linear medium), f0 = %.1f MHz\n', f0/1e6);
fprintf('  v (mm/s)  PRF (kHz)   delta (ns)   formula (dB)   simulated (dB)\n');
for vv = [1 10 50]
  for prf = [5e3]
    d = 2*vv*1e-3/c/prf; s2 = -interp1(t, s, t - d, 'linear', 0);
    fprintf('  %6.0f     %5.1f      %7.2f       %7.1f        %7.1f\n', vv, prf/1e3, d*1e9, ...
      20*log10(2*abs(sin(pi*f0*d))), dB10(Ef(s + s2)/Ef(s)));
  end
end

%% Band separation for filter-based harmonic imaging (Gaussian spectra)
% fundamental half-width at L dB: (B/2) f0 sqrt(L/6); the 2f0 spectrum of s^2 is sqrt(2) wider in Hz
for L = [6 20 30]
  fprintf('f0 and 2f0 spectra meet at -%2d dB when FBW = %.0f %%\n', L, 100*2/((1 + sqrt(2))*sqrt(L/6)));
end

%% Figure
figure('Position', [60 60 900 600]);
win = abs(t) < 1.2e-6;
subplot(2,2,1); plot(t(win)*1e6, eP(win), 'b', t(win)*1e6, eN(win), 'r'); grid on;
xlabel('t (\mus)'); ylabel('echo'); legend('+s', '-s'); title(sprintf('Tissue echoes, \\sigma = %.1f', sig));
subplot(2,2,2); plot(t(win)*1e6, pi_out(win), 'k', t(win)*1e6, am_out(win), 'm'); grid on;
xlabel('t (\mus)'); legend('PI sum', 'AM: e(s) - 2e(s/2)'); title('Combined outputs');
subplot(2,1,2); half = fk < 12e6; ref = max(abs(fft(eP)));
S1 = 20*log10(abs(fft(eP))/ref + 1e-9); S2 = 20*log10(abs(fft(pi_out))/ref + 1e-9); S3 = 20*log10(abs(fft(am_out))/ref + 1e-9);
plot(fk(half)/1e6, S1(half), 'b', fk(half)/1e6, S2(half), 'k', fk(half)/1e6, S3(half), 'm', 'LineWidth', 1.2);
ylim([-80 5]); grid on; xlabel('f (MHz)'); ylabel('dB re single-echo peak');
legend('single echo', 'PI sum', 'AM'); title('Spectra: PI removes f_0 and doubles 2f_0');
