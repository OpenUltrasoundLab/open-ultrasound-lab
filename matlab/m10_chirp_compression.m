%% m10_chirp_compression.m
% Linear FM (chirp) pulse compression: mainlobe width, range sidelobes,
% compression gain and the price of amplitude tapering.
%   - Rect-envelope LFM in complex baseband: the autocorrelation magnitude is
%     exactly (1-|tau|/T)|sinc(B tau (1-|tau|/T))|; the -6 dB width tends to
%     1.207/B and the peak sidelobe to -13.26 dB as TB grows.
%   - Compression gain at equal peak amplitude = energy ratio. Against a
%     Gaussian pulse with the same -6 dB bandwidth it is 1.505*TB.
%   - Tapers applied on receive only (mismatched filter), on transmit only,
%     or split as sqrt(w) on both sides.
% What to change: B, TB_list, TBw (taper table), the window list, Ntrial.
% Expected output: width*B 1.234 (TB 10) ... 1.208 (TB 200); PSL -14.9 ... -13.3 dB;
% Hann on receive at TB 50: PSL -31.3 dB, width*B 2.00, loss 1.76 dB; Hamming on
% receive: -36.5 dB (not -42.7: Fresnel ripple); Monte Carlo gain (10^4 trials) within about
% 0.25 dB of 10*log10(1.505*TB). Runtime about 3 s in Octave.
clear; close all; rng(1);

%% Helpers (anonymous only, so the script runs unchanged in MATLAB and Octave)
fs    = 200e6;                                    % complex-baseband sampling rate (Hz)
B     = 5e6;                                      % swept bandwidth (Hz)
tvec  = @(T) ((0:round(T*fs)-1)' - (round(T*fs)-1)/2)/fs;    % centred time axis (s)
lfm   = @(T) exp(1j*pi*(B/T)*tvec(T).^2);         % baseband up-chirp, instantaneous f = (B/T) t
sincf = @(x) (sin(pi*x) + (x==0))./(pi*x + (x==0));
nf    = @(n) 2^nextpow2(n);
% cross-correlation of received s with filter template h (output = s correlated with h), lag 0 centred
xcor  = @(s,h) fftshift(ifft(fft(s, nf(2*max(numel(s),numel(h)))).*conj(fft(h, nf(2*max(numel(s),numel(h)))))));
% -6 dB width (samples) with linear interpolation, around the peak index k0 of a normalized envelope a
xR    = @(a,i) i-1 + (a(i-1)-0.5)./(a(i-1)-a(i));
xL    = @(a,i) i+1 - (a(i+1)-0.5)./(a(i+1)-a(i));
w6f   = @(a,k0) xR(a, k0-1+find(a(k0:end) < 0.5, 1)) - xL(a, k0+1-find(a(k0:-1:1) < 0.5, 1));
% peak sidelobe (dB) outside the first minima on each side of the peak
pslf  = @(a,k0) 20*log10(max(a([1:(k0+1-find(diff(a(k0:-1:1)) > 0, 1))-1, (k0-1+find(diff(a(k0:end)) > 0, 1))+1:numel(a)])));
x01   = @(n) ((0:n-1)' + 0.5)/n;                  % window abscissa in (0,1)
wins  = {'rect','Tukey 25%','Hann','Hamming','Blackman'};
winf  = {@(x) ones(size(x)), ...
         @(x) 0.5 - 0.5*cos(2*pi*min(min(x,1-x),0.125)/0.25), ...   % Tukey: cosine ramps over 12.5% each side
         @(x) 0.5 - 0.5*cos(2*pi*x), ...
         @(x) 0.54 - 0.46*cos(2*pi*x), ...
         @(x) 0.42 - 0.5*cos(2*pi*x) + 0.08*cos(4*pi*x)};

%% 1. Rect LFM: width and PSL vs TB, against the closed form
TB_list = [10 20 50 100 200];
fprintf('Rect LFM, matched filter (complex baseband, B = %.0f MHz)\n   TB    width*B   PSL (dB)   closed-form max error\n', B/1e6);
for TBv = TB_list
    T = TBv/B; s = lfm(T);
    y = xcor(s,s); a = abs(y)/max(abs(y)); [~,k0] = max(a);
    tau = ((1:numel(a))' - k0)/fs;
    f = (abs(tau) < T).*(1 - abs(tau)/T).*abs(sincf(B*tau.*(1 - abs(tau)/T)));
    fprintf('%5d    %6.3f    %6.2f     %.1e\n', TBv, w6f(a,k0)/fs*B, pslf(a,k0), max(abs(a-f)));
end
fprintf('Limit TB -> inf: sinc, width*B = 1.207, PSL = -13.26 dB\n');

%% 2. Tapering: receive (mismatched), transmit, split; TB = 50
TBw = 50; T = TBw/B; s = lfm(T); n = numel(s); x = x01(n);
fprintf('\nTapers at TB = %d: PSL (dB) / width*B / SNR loss vs rect matched (dB)\n', TBw);
fprintf('  window       receive only            transmit only           split sqrt(w)\n');
ys = cell(1,numel(wins));
for k = 1:numel(wins)
    w = winf{k}(x); line = sprintf('  %-10s', wins{k});
    for mode = 1:3
        if mode == 1, st = s;           h = s.*w;        end     % mismatched receive filter
        if mode == 2, st = s.*w;        h = st;          end     % tapered transmit, matched receive
        if mode == 3, st = s.*sqrt(w);  h = st;          end     % split
        y = xcor(st,h); a = abs(y)/max(abs(y)); [~,k0] = max(a);
        loss = 10*log10(max(abs(y))^2/(sum(abs(h).^2)*n));     % rect matched reference: peak^2/sum|h|^2 = n
        line = [line sprintf('  %6.1f / %4.2f / %5.2f  ', pslf(a,k0), w6f(a,k0)/fs*B, loss)]; %#ok<AGROW>
        if mode == 1, ys{k} = a; end
    end
    fprintf('%s\n', line);
end
fprintf('Ideal spectral-window PSL for comparison: Hann -31.5, Hamming -42.7, Blackman -58.1 dB\n');

%% 3. Compression gain: Monte Carlo against a Gaussian pulse of the same bandwidth and peak amplitude
% A lower complex-baseband rate (5x the swept band) is enough here and keeps the Monte Carlo fast.
fs3 = 25e6;
lfm3 = @(T) exp(1j*pi*(B/T)*(((0:round(T*fs3)-1)' - (round(T*fs3)-1)/2)/fs3).^2);
sf = B/2.3548; st = 1/(2*pi*sf);                  % Gaussian amplitude spectrum std (Hz), envelope std (s)
tg = (-ceil(5*st*fs3):ceil(5*st*fs3))'/fs3; g = exp(-tg.^2/(2*st^2));   % reference pulse, peak 1
Ntrial = 10000; sigma = 3;                        % complex white noise, E|n|^2 = sigma^2
fprintf('\nCompression gain (equal peak amplitude 1, same -6 dB bandwidth)\n   TB   theory 1.505*TB   energy ratio   Monte Carlo\n');
for TBv = [10 50 100]
    T = TBv/B; s = lfm3(T);
    SNR = @(v) abs(mean(v))^2/mean(abs(v-mean(v)).^2);
    nz = sigma*(randn(numel(s),Ntrial) + 1j*randn(numel(s),Ntrial))/sqrt(2);
    yc = s'*(repmat(s,1,Ntrial) + nz);            % matched filter output at the peak
    nz = sigma*(randn(numel(g),Ntrial) + 1j*randn(numel(g),Ntrial))/sqrt(2);
    yg = g'*(repmat(g,1,Ntrial) + nz);
    fprintf('%5d     %6.2f dB        %6.2f dB      %6.2f dB\n', TBv, 10*log10(1.505*TBv), ...
        10*log10(sum(abs(s).^2)/sum(g.^2)), 10*log10(SNR(yc)/SNR(yg)));
end

%% Figures
T = 10e-6; s = lfm(T); n = numel(s); tt = tvec(T)*1e6;
figure('Position',[100 100 900 600]);
subplot(2,1,1);
f0 = 5e6; tr = tvec(T); plot(tr*1e6, real(s.*exp(1j*2*pi*f0*tr)), 'b'); xlim([-5 5]);
xlabel('time (\mus)'); ylabel('amplitude');
title(sprintf('RF chirp, f_0 = 5 MHz, B = %.0f MHz, T = %.0f \\mus (TB = %d)', B/1e6, T*1e6, round(T*B)));
subplot(2,1,2);
N = 2^16; S = fftshift(abs(fft(s, N))); fax = ((0:N-1)' - N/2)*fs/N;
plot(fax/1e6 + 5, S/max(S), 'b'); xlim([1 9]); grid on;
xlabel('frequency (MHz)'); ylabel('|S(f)| (normalized)');
title('Spectrum: nearly rectangular, Fresnel ripple at the band edges');

figure('Position',[100 100 900 420]);
lag = ((1:numel(ys{1}))' - (numel(ys{1})/2+1))/fs*1e6;
plot(lag, 20*log10(ys{1}+1e-12), 'k', lag, 20*log10(ys{3}+1e-12), 'b', lag, 20*log10(ys{4}+1e-12), 'r');
xlim([-10 10]); ylim([-80 3]); grid on;
xlabel('lag (\mus)'); ylabel('envelope (dB)');
legend('rect (matched)', 'Hann on receive', 'Hamming on receive', 'Location', 'northeast');
title(sprintf('Compressed chirp, TB = %d', TBw));
