%% m03_rf_filter_deconv.m
% Purpose : 1D pulse-echo line r = p * rho + noise at fs = 40 MHz. (1) Band-
%           pass filter it and show that the output is a sum of shifted,
%           scaled copies of the filtered pulse; compute the SNR gain for
%           three pass bands. (2) Wiener-deconvolve it and tabulate the
%           trade-off between effective pulse width and noise.
% Change  : sigma (noise std re pulse peak), the pass bands, the eps list.
% Expect  : copies-vs-filtered difference ~1e-16; expected SNR gains +5.8,
%           +7.3, +8.5 dB for 2-8, 3-7, 4-6 MHz; Wiener: max gain 1/(2 sqrt eps),
%           width 0.294 us without deconvolution, 0.214 us at eps = 1e-2,
%           0.155 us at 1e-4, with peak-SNR changes -0.6 and -15.7 dB.
clear; close all; rng(1);

fs = 40;                                    % sampling rate (MHz), time in us
f0 = 5; B = 3;                              % centre frequency, -6 dB bandwidth (MHz)
st = sqrt(2*log(2)) / (pi*B) * fs;          % Gaussian envelope std (samples)
k  = -20:20;                                % pulse support (samples)
p  = exp(-k.^2 / (2*st^2)) .* cos(2*pi*f0/fs*k);
N  = 256;                                   % 6.4 us line
pos = [50 90 96 150 200]; amp = [1 -0.7 0.6 0.5 0.8];
rho = zeros(1, N); rho(pos+1) = amp;
fprintf('pulse: sigma_t = %.3f samples, -6 dB envelope width %.3f us (%.3f mm)\n', st, ...
        2*sqrt(2*log(2))*st/fs, 2*sqrt(2*log(2))*st/fs*1.54/2);

sfull = conv(rho, p); s = sfull(21:20+N);   % centred: no delay
sigma = 0.1; w = sigma * randn(1, N); r = s + w;

%% Band-pass FIR (41-tap Hamming-windowed difference of sincs)
hw = 0.54 - 0.46*cos(2*pi*(0:40)/40);
sincd = @(f, m) (m == 0) * 2*f/fs + (m ~= 0) .* sin(2*pi*f/fs*m) ./ (pi*m + (m == 0));
bpf = @(f1, f2) hw .* (sincd(f2, k) - sincd(f1, k));
bands = [2 8; 3 7; 4 6];
for b = 1:3
    hb = bpf(bands(b,1), bands(b,2));
    yf = conv(r, hb);  y  = yf(21:20+N);    % centred output (20-sample delay removed)
    sf = conv(s, hb);  sf = sf(21:20+N);
    nf = conv(w, hb);  nf = nf(21:20+N);
    gS = sum(sf.^2) / sum(s.^2); gN = sum(hb.^2);
    fprintf('band %d-%d MHz: signal %+.2f dB, noise %+.2f dB, expected SNR gain %+.2f dB, this realization %+.2f dB\n', ...
            bands(b,1), bands(b,2), 10*log10(gS), 10*log10(gN), 10*log10(gS/gN), ...
            10*log10(sum(sf.^2)/sum(nf.^2)) - 10*log10(sum(s.^2)/sum(w.^2)));
    if b == 2, y37 = y; s37 = sf; hb37 = hb; end
end
q = conv(p, hb37);                           % filtered pulse, 81 taps, centre at index 41
copies = zeros(1, N);
for i = 1:numel(pos)
    idx = (pos(i)+1) + (-40:40);
    ok = idx >= 1 & idx <= N;
    copies(idx(ok)) = copies(idx(ok)) + amp(i) * q(ok);
end
fprintf('3-7 MHz: max |noise-free output - sum of shifted copies| = %.1e\n', max(abs(s37 - copies)));
fprintf('linear-phase delay of the 41-tap FIR: 20 samples = %.2f us = %.3f mm\n', 20/fs, 20/fs*1.54/2);

%% Wiener deconvolution: G = conj(P) / (|P|^2 + eps |P|max^2)
NF = 512;
pc = zeros(1, NF); pc(mod(k, NF) + 1) = p;    % pulse centred at sample 0 (circular)
P = fft(pc); Pm = max(abs(P));
fq = (0:NF-1) * fs / NF;
Pat = @(f) abs(sum(p .* exp(-1j*2*pi*f/fs*k)));   % DTFT at exactly f (MHz), not the nearest bin
fprintf('|P| re peak: %.1f dB at 2 MHz, %.1f dB at 8 MHz, %.1f dB at 10 MHz\n', ...
        20*log10(Pat(2)/Pat(5)), 20*log10(Pat(8)/Pat(5)), 20*log10(Pat(10)/Pat(5)));
R = fft(r, NF); Wn = fft(w, NF);
eps_list = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-7];
fprintf('%9s %10s %12s %12s %12s %14s\n', 'eps', 'max gain', 'width (us)', 'width (mm)', 'dSNR (dB)', 'noise std out');
est = zeros(numel(eps_list), N);
for i = 1:numel(eps_list)
    G = conj(P) ./ (abs(P).^2 + eps_list(i)*Pm^2);
    qe = real(ifft(P .* G));                  % effective pulse
    Q = fft(qe); Q(2:NF/2) = 2*Q(2:NF/2); Q(NF/2+2:end) = 0;
    e = fftshift(abs(ifft(Q))); e = e / max(e);
    c = find(e == 1, 1); il = c; ir = c;
    while e(il) >= 0.5, il = il - 1; end
    while e(ir) >= 0.5, ir = ir + 1; end
    wl = il + (0.5 - e(il)) / (e(il+1) - e(il));  wr = ir - (0.5 - e(ir)) / (e(ir-1) - e(ir));
    width = wr - wl;                           % samples
    rmsG = sqrt(mean(abs(G).^2));               % white-noise gain of G
    no = real(ifft(Wn .* G)); 
    xr = real(ifft(R .* G)); est(i, :) = xr(1:N);
    fprintf('%9.0e %+9.1f dB %12.3f %12.3f %+12.1f %14.3f\n', eps_list(i), 20*log10(max(abs(G))*Pm), ...
            width/fs, width/fs*1.54/2, 20*log10(qe(1)/rmsG), std(no(1:N)));
end

%% Figures
t = (0:N-1) / fs;
figure('Name', 'RF line: band-pass filtering');
subplot(3,1,1); stem(pos/fs, amp, 'k', 'filled'); xlim([0 N/fs]); ylim([-1.2 1.2]); title('reflectivity \rho');
subplot(3,1,2); plot(t, r, 'b'); xlim([0 N/fs]); title('RF line r = p * \rho + noise');
subplot(3,1,3); plot(t, copies, 'Color', [0.8 0.6 0.1]); hold on; plot(t, y37, 'Color', [0.1 0.5 0.4]); hold off;
xlim([0 N/fs]); xlabel('time (\mus)'); title('3-7 MHz output (green) and noise-free sum of copies (gold)');
figure('Name', 'Spectra and Wiener deconvolution');
subplot(2,1,1); half = 1:NF/2+1;
Hb = abs(fft([hb37(21:41) zeros(1, NF-41) hb37(1:20)]));
plot(fq(half), 20*log10(abs(P(half))/Pm), 'k', fq(half), 20*log10(Hb(half)), 'r', ...
     fq(half), 20*log10(abs(R(half))/max(abs(R))), 'b'); ylim([-80 5]); grid on;
xlabel('f (MHz)'); ylabel('dB'); title('|P| (k), |H_{bp}| (r), noisy |R| (b)');
subplot(2,1,2); plot(t, est(2,:), 'Color', [0.1 0.5 0.4]); hold on; plot(t, est(4,:) - 0.6, 'm');
stem(pos/fs, 0.3*amp, 'k'); hold off; xlim([0 N/fs]);
xlabel('time (\mus)'); title('Wiener: eps 1e-2 (green), 1e-4 (magenta, -0.6); 0.3 rho (k)');
