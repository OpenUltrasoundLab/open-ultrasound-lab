%% m04_fft_1d_basics.m  -- Module 04: frequency axes, aliasing, zero padding, Parseval, windows
% Purpose : build the FFT frequency axis correctly, see aliasing, show that zero padding
%           interpolates but does not resolve (a longer record does), check Parseval, and
%           measure the figures of merit of the classic windows.
% Change  : fs, N, the tone frequencies, the tone separation sepT (in units of 1/T),
%           the window list.
% Expect  : alias of 7 MHz at fs = 10 MHz -> 3.000 MHz; legal bandpass rates for HARD-supported
%           bands: [3.5, 6.5] MHz -> fs >= 13 or 6.5..7.0 MHz, [3, 7] MHz -> only fs >= 14 MHz
%           (band edges at 6.75 MHz land at 3.25 and 0.25 MHz); if [3.5, 6.5] are only the -6 dB
%           edges of a 5 MHz Gaussian, the zone-2 edges at 6.75 MHz sit at -7.07 / -8.19 dB, so
%           the tails alias; Hann tone vs noise normalization: -4.26 dB (noise, uncorrected) vs
%           -6.02 dB (tone), +1.76 dB left if CG^2 is used on noise; padded short record shows 1 peak,
%           4x longer record shows 2; Parseval gap ~1e-12; window table
%           rect -13.3 dB / Hann -31.5 dB / Hamming -42.7 dB / Blackman -58.1 dB,
%           roll-off about -6, -18, -6, -18 dB/octave.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

%% 1. Frequency axis and fftshift
fs = 40;                 % sampling rate (MHz)
N  = 1024;               % record length (samples)
T  = N/fs;               % record duration (us)
df = fs/N;               % bin spacing (MHz)
f_unshift = (0:N-1)*df;                   % bins as returned by fft (0 ... fs)
f_shift   = (-N/2:N/2-1)*df;              % axis after fftshift (-fs/2 ... fs/2)
n  = 0:N-1;
x  = cos(2*pi*5*n/fs + 0.3);              % 5 MHz tone
X  = fftshift(fft(x));
[~, ipk] = max(abs(X(N/2+1:end)));        % positive half only
fprintf('--- frequency axis ---\n');
fprintf('fs = %g MHz, N = %d: bin spacing fs/N = %.2f kHz, T = %.1f us, 1/T = %.2f kHz\n', fs, N, df*1e3, T, 1e3/T);
fprintf('5 MHz tone: peak at %.3f MHz (positive half of the fftshifted axis)\n', f_shift(N/2+ipk));
fprintf('unshifted bin k = N-1 is %.4f MHz = -df (negative frequency)\n', f_unshift(end) - fs);

%% 2. Aliasing: a 7 MHz tone sampled at 10 MHz
fsa = 10; Na = 4096; na = 0:Na-1;
w_hann = 0.5 - 0.5*cos(2*pi*na/Na);               % periodic Hann, written out (no toolbox)
xa = cos(2*pi*7*na/fsa + 0.6) .* w_hann;
Xa = abs(fft(xa)); [~, ka] = max(Xa(1:Na/2+1));
f_alias_formula = abs(7 - fsa*round(7/fsa));
fprintf('\n--- aliasing ---\n');
fprintf('7 MHz at fs = 10 MHz: formula |f0 - fs*round(f0/fs)| = %.3f MHz, FFT peak = %.3f MHz\n', ...
        f_alias_formula, (ka-1)*fsa/Na);

%% 2b. Bandpass sampling: for a band strictly confined to [fL, fH] (hard support),
%     2 fH/n <= fs <= 2 fL/(n-1), 1 <= n <= floor(fH/(fH-fL)); n = 1 has no upper bound
fprintf('\n--- bandpass sampling ---\n');
bands = [3.5 6.5; 3 7];                            % [fL fH] in MHz
for ib = 1:size(bands, 1)
    fL = bands(ib, 1); fH = bands(ib, 2); Bw = fH - fL; nmax = floor(fH/Bw);
    fprintf('band [%.1f, %.1f] MHz, B = %.1f, n_max = %d:', fL, fH, Bw, nmax);
    for nz = 1:nmax
        if nz == 1, fprintf('  n=1: fs >= %.2f', 2*fH);
        else, fprintf('  n=%d: %.2f..%.2f', nz, 2*fH/nz, 2*fL/(nz-1)); end
    end
    fprintf(' MHz\n');
end
fsb = 6.75; fe = [3.5 6.5];                        % band edges sampled at 6.75 MHz (zone n = 2)
fa_e = abs(fe - fsb*round(fe/fsb));
fprintf('fs = %.2f MHz: band edges %.1f / %.1f MHz alias to %.2f / %.2f MHz (zone 2: inverted)\n', fsb, fe, fa_e);
% A Gaussian has no hard support: if 3.5/6.5 MHz are its -6 dB edges, what lies outside zone 2?
sgf = 1.5/sqrt(2*log(2));                          % amplitude-spectrum std of a 5 MHz, 60 % Gaussian (MHz)
zed = [fsb/2 fsb];                                 % zone-2 edges at fs = 6.75 MHz
fprintf('5 MHz, 60 %% Gaussian: level at the zone-2 edges %.3f / %.2f MHz = %.2f / %.2f dB re peak\n', ...
        zed, 20*log10(exp(-(zed - 5).^2/(2*sgf^2))));
fprintf('(so -6 dB edges only mark where the -6 dB intervals do not overlap, not alias-free sampling)\n');

%% 3. Zero padding vs a longer record (complex tones, as IQ data)
f1 = 5.0; sepT = 0.6;                     % tone separation in units of 1/T
Ns = 64; Ts = Ns/fs; f2 = f1 + sepT/Ts;   % 0.6/T apart for the SHORT record
mk = @(Nr) exp(1j*2*pi*f1*(0:Nr-1)/fs) + exp(1j*2*pi*f2*(0:Nr-1)/fs);
countPeaks = @(m) sum( m(2:end-1) > m(1:end-2) & m(2:end-1) >= m(3:end) & m(2:end-1) > 0.5*max(m) );
Npad = 4096;
Xs_pad  = abs(fft(mk(Ns), Npad));         % short record, 64x zero padded
Xl_pad  = abs(fft(mk(4*Ns), Npad));       % 4x longer record, same display grid
fax = (0:Npad-1)*fs/Npad;
sel = fax > f1 - 1.5 & fax < f2 + 1.5;
fprintf('\n--- zero padding vs longer record ---\n');
fprintf('tones %.3f and %.3f MHz (%.2f/T apart for N = %d)\n', f1, f2, sepT, Ns);
fprintf('N = %3d padded to %d: %d peak(s) (bins of %.4f MHz, resolution 1/T = %.3f MHz)\n', Ns, Npad, countPeaks(Xs_pad(sel)), fs/Npad, 1/Ts);
fprintf('N = %3d padded to %d: %d peak(s) (resolution 1/T = %.3f MHz)\n', 4*Ns, Npad, countPeaks(Xl_pad(sel)), 1/(4*Ts));

%% 4. Parseval
xr = randn(1, N);
gap = sum(abs(xr).^2) - sum(abs(fft(xr)).^2)/N;
fprintf('\n--- Parseval ---\nsum|x|^2 - sum|X|^2/N = %.2e\n', gap);

%% 5. Window figures of merit (DFT-even windows, N = 256, 64x zero padding)
Nw = 256; P = 64; L = Nw*P; m = 0:Nw-1; a = 2*pi*m/Nw;
I0 = @(z) sum( bsxfun(@rdivide, bsxfun(@power, z(:)/2, 2*(0:40)), factorial(0:40).^2), 2 ).';  % Bessel I0 series
names = {'Rectangular', 'Hann', 'Hamming', 'Blackman', 'Kaiser b=6'};
W = [ ones(1, Nw);
      0.5 - 0.5*cos(a);
      0.54 - 0.46*cos(a);
      0.42 - 0.5*cos(a) + 0.08*cos(2*a);
      I0(6*sqrt(max(0, 1 - (2*m/Nw - 1).^2))) / I0(6) ];
theoryPSL = [-13.26 -31.47 -42.66 -58.11 NaN];
fprintf('\n--- window table (computed) ---\n');
fprintf('%-12s %8s %8s %8s %8s %6s %6s %8s %10s\n', 'window', 'PSL dB', '-3dB bin', '-6dB bin', 'null bin', 'ENBW', 'CG', 'scal dB', 'roll dB/oct');
figure('Name', 'windows', 'Position', [100 100 760 520]);
cols = lines(5);
for iw = 1:size(W, 1)
    w = W(iw, :);
    S = abs(fft(w, L)); S = S/S(1);
    d = 20*log10(S(1:L/2) + 1e-300);
    k = 2; while ~(d(k) < d(k-1) && d(k) <= d(k+1)), k = k + 1; end     % first null
    psl = max(d(k:end));
    w3 = 2*interp1(d(1:k), (0:k-1)/P, -3.0103);
    w6 = 2*interp1(d(1:k), (0:k-1)/P, -6.0206);
    enbw = Nw*sum(w.^2)/sum(w)^2;
    cg = mean(w);
    scal = -d(P/2 + 1);                                                % tone half a bin off
    env = @(b) max(d(b*P+1:(b+1)*P));                                  % sidelobe peak near bin b
    roll = (env(40) - env(10))/2;                                      % 2 octaves, per octave
    fprintf('%-12s %8.2f %8.3f %8.3f %8.3f %6.3f %6.3f %8.2f %10.1f', names{iw}, psl, w3, w6, (k-1)/P, enbw, cg, scal, roll);
    if ~isnan(theoryPSL(iw)), fprintf('   (theory PSL %.2f)\n', theoryPSL(iw)); else, fprintf('\n'); end
    plot((0:L/2-1)/P, d, 'Color', cols(iw, :), 'LineWidth', 1.1); hold on;
end
hold off; xlim([0 12]); ylim([-110 3]); grid on;
xlabel('frequency offset (bins of N)'); ylabel('|W| (dB re peak)');
title('Window spectra: mainlobe width vs sidelobe level'); legend(names, 'Location', 'northeast');

%% 5b. Tone vs noise normalization (periodic Hann): different corrections
wh = 0.5 - 0.5*cos(2*pi*(0:Nw-1)/Nw);
fprintf('\n--- Hann normalization (tone: divide by sum(w); PSD: |Xw|^2/(fs*sum(w.^2))) ---\n');
fprintf('CG = %.3f, mean(w^2) = %.4f, ENBW = %.3f bins\n', mean(wh), mean(wh.^2), Nw*sum(wh.^2)/sum(wh)^2);
fprintf('uncorrected tone amplitude: %.2f dB; uncorrected noise power: %.2f dB; noise corrected by CG^2: %+.2f dB\n', ...
        20*log10(mean(wh)), 10*log10(mean(wh.^2)), 10*log10(mean(wh.^2)/mean(wh)^2));

%% 6. Figure: aliasing and zero padding
figure('Name', 'aliasing and zero padding', 'Position', [100 100 760 600]);
subplot(2, 1, 1);
tc = linspace(0, 1.6, 2000); ts = 0:1/fsa:1.6;
plot(tc, cos(2*pi*7*tc + 0.6), 'Color', [0.6 0.6 0.6]); hold on;
plot(tc, cos(2*pi*3*tc - 0.6), 'r--', 'LineWidth', 1.4);
plot(ts, cos(2*pi*7*ts + 0.6), 'ko', 'MarkerFaceColor', 'k'); hold off;
xlim([0 1.6]); xlabel('time (us)'); ylabel('amplitude');
title('7 MHz tone (gray) sampled at 10 MHz (dots) = 3 MHz alias, phase -\phi (red)');
subplot(2, 1, 2);
plot(fax(sel), Xs_pad(sel)/max(Xs_pad), 'b', 'LineWidth', 1.3); hold on;
plot(fax(sel), Xl_pad(sel)/max(Xl_pad), 'r', 'LineWidth', 1.3);
plot([f1 f1], [0 1.05], 'k:', [f2 f2], [0 1.05], 'k:'); hold off;
xlabel('frequency (MHz)'); ylabel('|X| (norm.)');
title('Zero padding (blue, N = 64) vs 4x longer record (red, N = 256)');
legend('N=64, padded 64x', 'N=256, padded 16x', 'Location', 'northeast');
