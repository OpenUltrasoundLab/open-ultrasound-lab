%% m10_matched_filter_snr.m
% Matched filter output SNR in white Gaussian noise, three conventions, and
% the cost of using the wrong filter.
%   - Real RF samples, white noise of variance sigma^2: peak SNR = sum(s.^2)/sigma^2,
%     which is 2E/N0 for a two-sided noise PSD N0/2 (sigma^2 = (N0/2)*fs).
%   - Complex (analytic / IQ) data, SNR defined as |y|^2 / E|n_out|^2: 3 dB lower (E/N0).
%   - Gaussian band-pass filter of the wrong bandwidth (ratio r) or centre:
%     efficiency 2r/(1+r^2) and exp(-D^2/(sf^2+sh^2)) (narrowband theory).
% What to change: bw (fractional bandwidth), EsN0_dB, r_list, Ntrial.
% Expected output: measured RF SNR within 0.1 dB of theory; analytic |y|^2 SNR
% 3.0 dB lower; Re{} of the analytic output back at the RF value; r = 0.5 costs 0.97 dB.
% Runtime about 6 s in Octave.
clear; close all; rng(1);

%% Pulse and noise
fs = 40e6;              % sampling rate (Hz)
f0 = 5e6;               % centre frequency (Hz)
bw = 0.60;              % -6 dB fractional bandwidth of the amplitude spectrum
sf = bw*f0/2.3548;      % std of the Gaussian amplitude spectrum (Hz)
st = 1/(2*pi*sf);       % std of the Gaussian envelope (s)
t  = (-3e-6:1/fs:3e-6)';
s0 = exp(-t.^2/(2*st^2)).*cos(2*pi*f0*t);
EsN0_dB = 10;                            % target sum(s.^2)/sigma^2 in dB
sigma = 1;
s = s0*sqrt(10^(EsN0_dB/10)*sigma^2/sum(s0.^2));
SNRth = sum(s.^2)/sigma^2;               % = 2E/N0
fprintf('Pulse: %d samples, sum(s^2)/sigma^2 = %.2f dB (theory 2E/N0)\n', numel(s), 10*log10(SNRth));
fprintf('Gain over the best single sample: %.2f dB\n', 10*log10(sum(s.^2)/max(s.^2)));

%% Monte Carlo: matched filter on RF and on the analytic signal (and mismatched filters)
Nchunk = 2; Nper = 4000; Ntrial = Nchunk*Nper;  % 8000 trials in chunks to save memory
Nrec = 512; i0 = 256;                          % record length, echo centre index
x0 = zeros(Nrec,1); L = numel(s); x0(i0-(L-1)/2 : i0+(L-1)/2) = s;
win = i0-(L-1)/2 : i0+(L-1)/2;
wA = zeros(Nrec,1); wA(1) = 1; wA(2:Nrec/2) = 2; wA(Nrec/2+1) = 1;   % analytic-signal weights
sA = ifft(fft(x0).*wA); sA = sA(win);
r_list = [0.25 0.5 0.75 1 1.5 2 3 4];
effTh = 2*r_list./(1+r_list.^2); effEx = zeros(size(r_list));
Hs = cell(size(r_list)); rows = zeros(size(r_list)); Nf = Nrec + 512;
for k = 1:numel(r_list)
    sth = st/r_list(k);                          % filter envelope std
    th = (-ceil(6*sth*fs):ceil(6*sth*fs))'/fs;
    h = exp(-th.^2/(2*sth^2)).*cos(2*pi*f0*th);  % symmetric, so h(-t) = h(t)
    ys = conv(s, h); [~, kp] = max(abs(ys));
    effEx(k) = ys(kp)^2/(sigma^2*sum(h.^2))/SNRth;   % exact deterministic efficiency
    Hs{k} = fft(h, Nf); rows(k) = win(1) + kp - 1;  % output index of the signal peak
end
yRF = zeros(1,Ntrial); yA = zeros(1,Ntrial); yK = zeros(numel(r_list),Ntrial);
for c = 1:Nchunk
    X = repmat(x0,1,Nper) + sigma*randn(Nrec,Nper);
    idx = (c-1)*Nper + (1:Nper);
    yRF(idx) = s' * X(win,:);                    % matched filter output at the peak
    XA = ifft(fft(X).*repmat(wA,1,Nper));
    yA(idx) = sA' * XA(win,:);                   % sA' is the conjugate transpose
    Xf = fft([X; zeros(Nf-Nrec,Nper)]);
    for k = 1:numel(r_list)
        Y = real(ifft(Xf.*repmat(Hs{k},1,Nper)));
        yK(k,idx) = Y(rows(k),:);
    end
end
SNRrf = mean(yRF)^2/var(yRF);
mA = mean(yA);
SNRc  = abs(mA)^2/mean(abs(yA-mA).^2);
SNRre = real(mA)^2/var(real(yA));
effMC = (mean(yK,2).^2./var(yK,0,2))'/SNRth;
fprintf('\nMatched filter, %d trials\n', Ntrial);
fprintf('  RF (real), coherent        : %.2f dB  (theory %.2f)\n', 10*log10(SNRrf), 10*log10(SNRth));
fprintf('  analytic, |y|^2/E|n|^2     : %.2f dB  (theory %.2f = E/N0)\n', 10*log10(SNRc), 10*log10(SNRth/2));
fprintf('  analytic, Re{y} (coherent) : %.2f dB\n', 10*log10(SNRre));

fprintf('\nGaussian filter of bandwidth ratio r (filter/echo):\n   r    theory 2r/(1+r^2)   exact     Monte Carlo\n');
for k = 1:numel(r_list)
    fprintf('%5.2f   %7.2f dB        %7.2f dB  %7.2f dB\n', r_list(k), 10*log10(effTh(k)), 10*log10(effEx(k)), 10*log10(effMC(k)));
end
for D = [0.5 1 2]
    h = exp(-t.^2/(2*st^2)).*cos(2*pi*(f0+D*sf)*t);
    ys = conv(s,h); Nz = 2048; wz = [1; 2*ones(Nz/2-1,1); 1; zeros(Nz/2-1,1)];
    A = abs(ifft(fft(ys, Nz).*wz));               % envelope of the filter output
    fprintf('centre offset %.1f sf (%.2f MHz): loss %.2f dB, theory exp(-D^2/2) = %.2f dB\n', D, D*sf/1e6, ...
        10*log10(max(A)^2/(sigma^2*sum(h.^2))/SNRth), 10*log10(exp(-D^2/2)));
end

%% Figures
figure('Position',[100 100 900 560]);
subplot(2,1,1);
tt = (0:Nrec-1)'/fs*1e6; X = repmat(x0,1,1) + sigma*randn(Nrec,1);
plot(tt, X(:,1), 'Color', [0.6 0.6 0.6]); hold on; plot(tt, x0, 'b', 'LineWidth', 1.5);
xlim([2 11]); xlabel('time (\mus)'); ylabel('amplitude'); title('Echo in white noise (one trial)');
subplot(2,1,2);
y1 = conv(X(:,1), flipud(s), 'same');
y0 = conv(x0, flipud(s), 'same');                % noise-free output
plot(tt, y1/sqrt(sigma^2*sum(s.^2)), 'k'); hold on;
plot(tt, y0/sqrt(sigma^2*sum(s.^2)), 'b', 'LineWidth', 1.2);
plot([0 tt(end)], sqrt(SNRth)*[1 1], 'r--');
legend('signal + noise', 'signal only', 'sqrt(2E/N0)', 'Location', 'southeast');
xlim([2 11]); xlabel('time (\mus)'); ylabel('output / noise std');
title(sprintf('Matched filter output: peak SNR theory %.1f dB', 10*log10(SNRth)));

figure('Position',[100 100 700 420]);
semilogx(r_list, 10*log10(effTh), 'b-', r_list, 10*log10(effEx), 'ko', r_list, 10*log10(effMC), 'r+', 'LineWidth', 1.2);
grid on; xlabel('filter bandwidth / echo bandwidth  r'); ylabel('SNR loss (dB)');
legend('2r/(1+r^2) narrowband theory', 'exact', 'Monte Carlo', 'Location', 'south');
title('Matched filters are flat near the optimum');
