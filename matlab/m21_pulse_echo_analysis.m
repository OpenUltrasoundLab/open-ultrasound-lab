%% m21_pulse_echo_analysis.m
% Pulse-echo test analysis on a synthetic "measured" echo from a flat steel reflector.
% What it shows
%   centre frequency (arithmetic and geometric), -6 dB bandwidth, -6 / -20 dB pulse length,
%   loop sensitivity and its correction for the reflector and water path, and how the gate
%   length changes the answer when a lens reverberation is present.
% What to change: fc, bw (two-way -6 dB FBW), reverberation level/delay, SNR, gate length.
% Expected output (defaults): fc ~ 5.0 MHz, FBW ~ 65 %, loop sensitivity ~ -20.9 dB for a true
%   two-way insertion loss of -20 dB (steel -0.57 dB, water -0.33 dB); FBW jumps once the gate
%   includes the reverberation (-14 dB at 1.6 us) and the spectrum ripples with period 1/1.6 us.
% Units: us, MHz, V.  Base MATLAB / Octave.

clear; close all; rng(1);
fs = 250;  N = 2000;  t = (0:N-1)/fs;   % 250 MS/s scope, 8 us record
fc = 5;  bw = 0.65;                     % transducer: centre frequency (MHz), two-way -6 dB FBW
rv_dB = -14;  dl = 1.6;                 % lens reverberation level (dB) and delay (us)
snr_dB = 40;  IL = -20;                 % echo peak / rms noise; true two-way insertion loss (dB)
Vtx = 100;  z = 30;                     % excitation (Vpp), reflector distance (mm)
Rst = (45.4 - 1.48)/(45.4 + 1.48);      % steel / water amplitude reflection
aw  = 2.2e-4;                           % water attenuation, dB/(mm MHz^2)
TE  = 1.0;                              % echo arrival within the record (us)

% causal gamma-tone echo, n = 3: -6 dB FBW = 0.2049/(fc*tau)
n = 3; tau = 0.2049/(bw*fc);
g = @(tt) (tt > 0) .* (max(tt,0)/tau).^n .* exp(-max(tt,0)/tau) .* cos(2*pi*fc*tt);
e0 = g(t - TE);
Vpp = Vtx * 10^(IL/20) * Rst * 10^(-2*z*aw*fc^2/20);
x = (e0 - 10^(rv_dB/20) * g(t - TE - dl)) * Vpp / (max(e0) - min(e0));
x = x + (Vpp/2) * 10^(-snr_dB/20) * randn(size(x));

% analysis function (anonymous helpers only)
nfft = 16384;
tukey = @(L, nt) [0.5-0.5*cos(pi*(0:nt-1)/nt), ones(1, L-2*nt), 0.5-0.5*cos(pi*(nt-1:-1:0)/nt)];

gates = [1.2 0.6 0.8 1.0 1.6 2.0 2.5 3.0];
fprintf('gate (us)   fpk    fl     fu     fc(ar)  fc(geo)  FBW(%%)  d-6(us)  d-20(us)\n');
for ig = 1:numel(gates)
    gl = gates(ig);
    i0 = round((TE - 0.25)*fs) + 1; L = round(gl*fs); nt = max(2, round(min(0.08*L, 0.1*fs)));
    seg = x(i0:i0+L-1) .* tukey(L, nt);
    X = fft(seg, nfft); mag = abs(X(1:nfft/2+1)); df = fs/nfft; f = (0:nfft/2)*df;
    [mx, ip] = max(mag(2:end)); ip = ip + 1; thr = mx/2;
    kl = ip; while kl > 2 && mag(kl) >= thr, kl = kl - 1; end
    ku = ip; while ku < numel(mag) && mag(ku) >= thr, ku = ku + 1; end
    fl = f(kl) + (thr - mag(kl))/(mag(kl+1) - mag(kl))*df;
    fu = f(ku-1) + (mag(ku-1) - thr)/(mag(ku-1) - mag(ku))*df;
    fcA = (fl + fu)/2; fcG = sqrt(fl*fu); fbw = (fu - fl)/fcA;
    % envelope by the FFT analytic signal
    Xa = X; Xa(2:nfft/2) = 2*Xa(2:nfft/2); Xa(nfft/2+2:end) = 0; xa = ifft(Xa); env = abs(xa(1:L));
    [em, ie] = max(env);
    lev = @(dbv) em*10^(dbv/20);
    ia = ie; while ia > 1 && env(ia) >= lev(-6), ia = ia - 1; end
    ib = ie; while ib < L && env(ib) >= lev(-6), ib = ib + 1; end
    d6 = (ib - ia)/fs;
    ia = ie; while ia > 1 && env(ia) >= lev(-20), ia = ia - 1; end
    ib = ie; while ib < L && env(ib) >= lev(-20), ib = ib + 1; end
    d20 = (ib - ia)/fs;
    fprintf('%6.2f    %5.2f  %5.2f  %5.2f   %5.2f   %5.2f    %5.1f    %5.3f    %5.3f\n', gl, f(ip), fl, fu, fcA, fcG, 100*fbw, d6, d20);
    if ig == 1
        LS = 20*log10((max(seg) - min(seg))/Vtx);
        ILc = LS - 20*log10(Rst) + 2*z*aw*fcA^2;
        d6P = d6; segP = seg; envP = env; tP = t(i0:i0+L-1); fP = f; magP = mag; fl1 = fl; fu1 = fu; fc1 = fcA; fbw1 = fbw;
    end
end
fprintf('\nGaussian time-bandwidth check: 0.883/(FBW*fc) = %.3f us; the causal gamma-tone echo measures %.3f us\n', 0.8825/(fbw1*fc1), d6P);
fprintf('Loop sensitivity (gate 1.2 us) = %.2f dB re %d Vpp\n', LS, Vtx);
fprintf('Corrected for steel (%.2f dB) and water (%.2f dB): %.2f dB (true %.1f dB)\n', 20*log10(Rst), -2*z*aw*fc1^2, ILc, IL);
fprintf('Arithmetic - geometric centre = %.3f MHz; fc*FBW^2/8 = %.3f MHz\n', fc1 - sqrt(fl1*fu1), fc1*fbw1^2/8);

figure(1); clf;
subplot(2,1,1); plot(t, x, 'b'); hold on; plot(tP, envP, 'Color', [0 0.5 0.4], 'LineWidth', 1.5);
yl = [-0.7 0.7]*Vpp; plot([tP(1) tP(1)], yl, '--', 'Color', [0.8 0.55 0]); plot([tP(end) tP(end)], yl, '--', 'Color', [0.8 0.55 0]); ylim(yl);
xlabel('t (us)'); ylabel('V'); title('Measured-like echo with gate (1.2 us) and envelope'); xlim([0 5]); grid on;
subplot(2,1,2); plot(fP, 20*log10(magP/max(magP)), 'b', 'LineWidth', 1.5); hold on;
plot([0 2*fc], [-6.02 -6.02], 'r--'); plot([fl1 fu1], [-6.02 -6.02], 'ro', 'MarkerFaceColor', 'r');
xlim([0 2*fc]); ylim([-40 2]); xlabel('f (MHz)'); ylabel('dB re peak'); grid on;
title(sprintf('Spectrum: fc = %.2f MHz, FBW = %.1f %%', fc1, 100*fbw1));
