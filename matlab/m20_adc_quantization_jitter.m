%% m20_adc_quantization_jitter.m
% ADC quantization SNR, aperture-jitter SNR, ENOB and processing gain (Module 20.5).
%   A coherently sampled sine (integer, odd number of cycles in the record)
%   is quantized to N bits and/or sampled with Gaussian clock jitter sigma_t.
%   SNR is measured from the FFT as signal-bin power over all other bins (SINAD).
% What to change: bit counts, jitter values, input frequency, input level.
% Expected output: measured quantization SNR within about 0.3 dB of
%   6.02N + 1.76; jitter-only SNR at 10 MHz near 84.0 dB (1 ps) and 70.1 dB
%   (5 ps); 14 bit + 1 ps at 10 MHz about 82 dB (ENOB about 13.3);
%   in-band SNR in a 4 MHz band at 40 MS/s about 7 dB above the Nyquist-band value.
clear; close all; rng(1);

%% Setup
fs  = 65e6;                 % sample rate (Hz)
Ns  = 8192;                 % record length
n   = 0:Ns-1;
k   = round(10e6/fs*Ns); if mod(k,2) == 0, k = k + 1; end
fin = k*fs/Ns;              % coherent input frequency (Hz), about 10 MHz
A   = 10^(-0.1/20);         % input level: -0.1 dBFS (full scale = +-1)
quant = @(x,N) max(-1, min(1 - 2^(1-N), round(x*2^(N-1))/2^(N-1)));
sinad = @(X,kk) 10*log10( abs(X(kk+1))^2 / (sum(abs(X(2:Ns/2)).^2) - abs(X(kk+1))^2) );   % X = fft(x)
ph0 = 0.37;
fprintf('fin = %.4f MHz (bin %d of %d), fs = %.1f MS/s, level %.1f dBFS\n', fin/1e6, k, Ns, fs/1e6, 20*log10(A));

%% Quantization only
fprintf('\n%5s %12s %12s %8s\n', 'bits', 'measured dB', 'theory dB', 'ENOB');
bits = [8 10 12 14 16]; mq = zeros(size(bits));
for j = 1:numel(bits)
  x = quant(A*sin(2*pi*fin*n/fs + ph0), bits(j));
  mq(j) = sinad(fft(x), k);
  th = 6.02*bits(j) + 1.76 + 20*log10(A);
  fprintf('%5d %12.2f %12.2f %8.2f\n', bits(j), mq(j), th, (mq(j) - 20*log10(A) - 1.76)/6.02);
end

%% Jitter only (no quantization)
fprintf('\n%10s %12s %12s\n', 'jitter ps', 'measured dB', 'theory dB');
sj = [0.5 1 2 5]*1e-12; mj = zeros(size(sj));
for j = 1:numel(sj)
  x = A*sin(2*pi*fin*(n/fs + sj(j)*randn(1,Ns)) + ph0);
  mj(j) = sinad(fft(x), k);
  fprintf('%10.1f %12.2f %12.2f\n', sj(j)*1e12, mj(j), -20*log10(2*pi*fin*sj(j)));
end

%% Combined: 14 bits + 1 ps at about 10 MHz
x = quant(A*sin(2*pi*fin*(n/fs + 1e-12*randn(1,Ns)) + ph0), 14);
mc = sinad(fft(x), k);
thc = -10*log10(10^(-(86.04 + 20*log10(A))/10) + 10^(20*log10(2*pi*fin*1e-12)/10));
fprintf('\n14 bit + 1 ps at %.2f MHz: measured %.2f dB, theory %.2f dB, ENOB %.2f\n', fin/1e6, mc, thc, (mc - 1.76)/6.02);
fprintf('jitter needed for 86 dB at 15 MHz: %.2f ps\n', 10^(-86/20)/(2*pi*15e6)*1e12);

%% Processing gain: 12 bits at 40 MS/s, 5 MHz tone, 4 MHz band
fs2 = 40e6; k2 = round(5e6/fs2*Ns); if mod(k2,2) == 0, k2 = k2 + 1; end
x2  = quant(A*sin(2*pi*k2/Ns*n + ph0), 12);
X2  = abs(fft(x2)).^2; fb = (0:Ns-1)*fs2/Ns;
band = find(fb >= 3e6 & fb <= 7e6); band = band(band ~= k2+1);
snrN = sinad(fft(x2), k2);
snrB = 10*log10(X2(k2+1)/sum(X2(band)));
fprintf('\n12 bit, 40 MS/s: SNR over 0-20 MHz %.2f dB; in 3-7 MHz %.2f dB; gain %.2f dB (theory 10log10(fs/2B) = %.2f dB)\n', ...
  snrN, snrB, snrB - snrN, 10*log10(fs2/2/4e6));
fprintf('LSB (2 Vpp): 12 bit %.1f uV, 14 bit %.1f uV\n', 2/2^12*1e6, 2/2^14*1e6);

%% Figures
figure('Name','m23 ADC spectra', 'Position', [50 50 760 720]);
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 7.5]);
x = quant(A*sin(2*pi*fin*(n/fs + 1e-12*randn(1,Ns)) + ph0), 12);
P = abs(fft(x)).^2*4/Ns^2; fa = (0:Ns/2-1)*fs/Ns/1e6;
subplot(2,1,1); plot(fa, 10*log10(P(1:Ns/2) + 1e-30)); grid on; ylim([-140 5]);
xlabel('frequency (MHz)'); ylabel('dBFS per bin'); title('12 bit, 1 ps jitter, 10 MHz tone');
subplot(2,1,2);
ff = linspace(0.5e6, 30e6, 200);
plot(ff/1e6, -20*log10(2*pi*ff*1e-12), ff/1e6, -20*log10(2*pi*ff*5e-12), ff/1e6, 74 + 0*ff, '--', ff/1e6, 86.04 + 0*ff, '--', fin/1e6, mj(2), 'o', fin/1e6, mj(4), 's');
grid on; ylim([50 110]); xlabel('input frequency (MHz)'); ylabel('SNR (dB)');
legend('1 ps', '5 ps', '12 bit', '14 bit', 'measured 1 ps', 'measured 5 ps', 'Location', 'southwest'); title('jitter-limited SNR');
figure('Name','m23 quantization SNR');
plot(bits, mq, 'o', 6:16, 6.02*(6:16) + 1.76, '-'); grid on;
xlabel('bits N'); ylabel('SNR (dB)'); legend('measured', '6.02N + 1.76', 'Location', 'northwest'); title('quantization SNR vs bits');
