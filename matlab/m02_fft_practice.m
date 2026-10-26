%% m02_fft_practice.m  -- Module 02: FFTs in practice
% Purpose : build frequency axes for fft and fftshift output (even and odd N), read a one-sided
%           amplitude spectrum in MHz, compute two real FFTs with one complex FFT, verify that
%           fft2 is fft over rows then columns, and measure the coherent processing gain N/2.
% Change  : fs, f0, N, the noise level sigma, the tone amplitude A.
% Expect  : 5 MHz pulse peaks at 5.00 MHz on both axes; fftshift puts DC at index floor(N/2)+1;
%           ifftshift(fftshift(X)) == X for odd N while fftshift(fftshift(X)) does not;
%           two-for-one and rows-then-columns errors ~1e-14; bin SNR within a few tenths of a dB
%           of 10log10(N/2) + per-sample SNR (A = 0.5, sigma = 1: -9.03 + 27.09 = 18.06 dB at N = 1024).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

%% 1. Frequency axes for an RF pulse
fs = 40;                         % sampling rate (MHz)
f0 = 5;                          % centre frequency (MHz)
N  = 1024;                       % FFT length
t  = ((0:N-1) - N/2)/fs;         % time (us), pulse centred in the record
x  = exp(-(t/0.25).^2) .* cos(2*pi*f0*t);   % Gaussian-envelope RF pulse, 1/e half-width 0.25 us
X  = fft(x);
f_fft   = (0:N-1)*fs/N;                          % fft order
f_shift = (-floor(N/2):ceil(N/2)-1)*fs/N;        % after fftshift
Xs = fftshift(X);
[~, i1] = max(abs(X(1:N/2+1)));  [~, i2] = max(abs(Xs) .* (f_shift > 0));
fprintf('--- frequency axes (N = %d, fs = %g MHz, df = %.2f kHz) ---\n', N, fs, 1e3*fs/N);
fprintf('peak on fft axis: %.2f MHz (index %d); on shifted axis: %.2f MHz\n', f_fft(i1), i1, f_shift(i2));
fprintf('DC after fftshift sits at index %d = floor(N/2)+1\n', find(f_shift == 0));
A1 = abs(X(1:N/2+1))/N; A1(2:end-1) = 2*A1(2:end-1);   % one-sided amplitude
fprintf('one-sided amplitude peak = %.4f (pulse spectrum, not a tone)\n', max(A1));

%% 2. Odd N: fftshift vs ifftshift
No = 7; v = 0:No-1;
fprintf('\n--- odd N = %d ---\n', No);
fprintf('fftshift(0:6)            = [%s]\n', num2str(fftshift(v)));
fprintf('centred axis (bins)      = [%s]\n', num2str(-floor(No/2):ceil(No/2)-1));
fprintf('ifftshift(fftshift(v))   = [%s]  (restores)\n', num2str(ifftshift(fftshift(v))));
fprintf('fftshift(fftshift(v))    = [%s]  (does not)\n', num2str(fftshift(fftshift(v))));

%% 3. Two real FFTs with one complex FFT
a = randn(1, N); b = randn(1, N);
Z = fft(a + 1j*b);
Zc = conj(Z([1, N:-1:2]));                       % Z*[N-k], indices mod N
Xa = (Z + Zc)/2; Xb = (Z - Zc)/(2j);
fprintf('\n--- two-for-one ---\nmax |Xa - fft(a)| = %.1e, max |Xb - fft(b)| = %.1e\n', max(abs(Xa - fft(a))), max(abs(Xb - fft(b))));

%% 4. fft2 = rows then columns
img = randn(64, 128);
R = fft(img, [], 2);                             % every row
C = fft(R, [], 1);                               % then every column
fprintf('max |fft2 - rows-then-columns| = %.1e\n', max(max(abs(fft2(img) - C))));

%% 5. Coherent processing gain
A = 0.5; sigma = 1; trials = 400;           % tone at fs/8, i.e. bin N/8 (always on a bin)
snr_sample = 10*log10(A^2/2/sigma^2);
fprintf('\n--- coherent gain (A = %.1f, sigma = %.0f, per-sample SNR %.2f dB) ---\n', A, sigma, snr_sample);
fprintf('     N   theory bin SNR (dB)   simulated (dB)\n');
for Nn = [64 256 1024 4096]
    nn = 0:Nn-1; ks = Nn/8 + 1;                          % MATLAB index of bin N/8
    s = A*cos(2*pi*(Nn/8)*nn/Nn);
    P_noise = 0;
    for r = 1:trials
        Wn = fft(sigma*randn(1, Nn));
        P_noise = P_noise + mean(abs(Wn([2:Nn/2])).^2);   % noise power per bin (positive half, no DC)
    end
    P_noise = P_noise/trials;
    Ssig = abs(fft(s)); P_tone = Ssig(ks)^2;
    fprintf('%6d   %18.2f   %14.2f\n', Nn, snr_sample + 10*log10(Nn/2), 10*log10(P_tone/P_noise));
end

%% 6. Figure
figure('Position', [100 100 900 600]);
subplot(2, 1, 1); plot(f_fft, abs(X)); xlim([0 fs]); grid on;
xlabel('frequency (MHz), fft order'); ylabel('|X|'); title('fft order: positive half, then the mirrored negative frequencies');
subplot(2, 1, 2); plot(f_shift, abs(Xs)); xlim([-fs/2 fs/2]); grid on;
xlabel('frequency (MHz), after fftshift'); ylabel('|X|'); title('fftshift: DC in the middle, peaks at \pm5 MHz');
