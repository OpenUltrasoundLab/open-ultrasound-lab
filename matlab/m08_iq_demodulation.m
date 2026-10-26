%% m08_iq_demodulation.m  -  IQ demodulation: mix, low-pass, decimate
% Purpose : build the IQ chain by hand: complex mixing with exp(-j2*pi*fd*t),
%           a windowed-sinc FIR low-pass you design yourself (Hamming window
%           written out, no toolbox), and decimation by M. Shows the spectrum
%           at each stage, compares 2|IQ| with the analytic-signal envelope,
%           sweeps M to show where aliasing starts, and shows what a
%           mixing-frequency mismatch does (residual rotation, band clipping).
% Change  : f0, bw, fd (demodulation frequency), fc (LPF cutoff), L (taps), M.
% Expect  : M = 4 (10 MHz IQ): envelope error ~0.6 % (set by the LPF);
%           reconstruction error stays at that floor to M = 5 (8 MHz) and
%           grows from M = 6; ~4 % at M = 8 (5 MHz), ~12 % at M = 9-10.
%           fd = f0 + 1 MHz rotates IQ at -1 MHz with the envelope intact
%           (1.7 %) while fc = 4 MHz passes the shifted band; fc = 2 MHz
%           clips it (envelope error ~24 %, rotation estimate biased).
clear; close all; rng(1);

%% Parameters
fs  = 40;  N = 2048;  t = (0:N-1)'/fs;     % MHz, samples, us
f0  = 5;   bw = 0.6;                       % echo center (MHz), -6 dB frac. BW
sigf = bw*f0/(2*sqrt(2*log(2)));  sigt = 1/(2*pi*sigf);
fd  = 5;                                   % demodulation frequency (MHz)
fc  = 4;                                   % LPF cutoff (MHz)
L   = 49;                                  % FIR length (odd -> integer delay)
M   = 4;                                   % decimation factor

% RF line: a few echoes of different phase plus weak noise
tc  = [8 11 11.4 14 20];  amp = [1 0.6 -0.8 0.4 0.7];
x   = zeros(N,1);
for k = 1:numel(tc)
  x = x + amp(k)*exp(-(t-tc(k)).^2/(2*sigt^2)).*cos(2*pi*f0*(t-tc(k)) + k);
end
x = x + 0.002*randn(N,1);

% Reference: analytic signal by FFT mask
h  = zeros(N,1); h(1) = 1; h(2:N/2) = 2; h(N/2+1) = 1;
xa = ifft(fft(x).*h);

% Windowed-sinc low-pass (Hamming window written out), unit DC gain
% sinc without the toolbox: sin(pi x)/(pi x), with sinc(0) = 1
sinc_ = @(u) (sin(pi*u) + (u==0)) ./ (pi*u + (u==0));
designLPF = @(fc_, L_) ( (2*fc_/fs)*sinc_((-(L_-1)/2:(L_-1)/2)'*2*fc_/fs) ...
                        .* (0.54 - 0.46*cos(2*pi*(0:L_-1)'/(L_-1))) );
normDC   = @(b) b/sum(b);
lpf      = @(v, b) conv(v, b, 'same');     % zero-phase alignment for odd L

% Full chain as one anonymous function: returns IQ at fs (before decimation)
iqchain  = @(x_, fd_, fc_) 2*lpf(x_.*exp(-1j*2*pi*fd_*t), normDC(designLPF(fc_, L)));

%% 1. The chain at nominal settings
b    = normDC(designLPF(fc, L));
xm   = x.*exp(-1j*2*pi*fd*t);              % stage 2: mixed
xl   = 2*lpf(xm, b);                       % stage 3: filtered (x2 restores amplitude)
iq   = xl(1:M:end);  fsIQ = fs/M;          % stage 4: decimated
ref  = xa.*exp(-1j*2*pi*fd*t);             % ideal IQ = analytic signal, carrier removed

% Envelope error at the kept samples
idx  = (1:M:N)';  sel = abs(xa(idx)) > 0.05;
eKept = 100*sqrt(sum((abs(iq(sel)) - abs(xa(idx(sel)))).^2)/sum(abs(xa(idx(sel))).^2));
% Band-limited reconstruction back to 40 MHz (what an interpolator sees)
Nd   = numel(iq);  rec = interpft(iq, Nd*M);
eRec = 100*sqrt(sum(abs(abs(rec(1:N)) - abs(xa)).^2.*(abs(xa)>0.05))/sum(abs(xa).^2.*(abs(xa)>0.05)));
% Complex error vs ideal IQ (includes phase)
ws    = abs(xa) > 0.05;                 % evaluate where the echoes are (noise floor excluded)
eCplx = 100*norm(xl(ws) - ref(ws))/norm(ref(ws));
fprintf('fd = %.1f MHz, fc = %.1f MHz, L = %d, M = %d (IQ rate %.1f MHz)\n', fd, fc, L, M, fsIQ);
fprintf('  2|IQ| vs |x_a| at kept samples      = %.3f %% rms\n', eKept);
fprintf('  band-limited reconstruction error   = %.3f %% rms\n', eRec);
fprintf('  complex IQ error before decimation  = %.3f %% rms\n', eCplx);
fprintf('  -20 dB half-width of the echo band  = %.2f MHz (needs fs/M > %.1f MHz)\n', ...
        sigf*sqrt(2*log(10)), 2*sigf*sqrt(2*log(10)));

%% 2. Sweep the decimation factor
Ms = 1:12; eM = zeros(size(Ms));
for i = 1:numel(Ms)
  q  = xl(1:Ms(i):end);  nd = numel(q);
  r  = interpft(q, nd*Ms(i));  r = r(1:min(N, nd*Ms(i)));
  rr = xa(1:numel(r));  w = abs(rr) > 0.05;
  eM(i) = 100*sqrt(sum((abs(r(w)) - abs(rr(w))).^2)/sum(abs(rr(w)).^2));
end
fprintf('\nM :      '); fprintf('%6d', Ms); fprintf('\nfsIQ:    '); fprintf('%6.1f', fs./Ms);
fprintf('\nerr (%%): '); fprintf('%6.2f', eM); fprintf('\n');

%% 3. Mixing-frequency mismatch
fdm = 6;                                    % demodulate 1 MHz above f0
for fcm = [4 2]
  xlm  = iqchain(x, fdm, fcm);
  ephm = 100*sqrt(sum((abs(xlm) - abs(xa)).^2.*(abs(xa)>0.05))/sum(abs(xa).^2.*(abs(xa)>0.05)));
  % residual rotation rate from the lag-one product inside the first echo
  w8   = abs(t - tc(1)) < 2*sigt;
  rot  = angle(sum(conj(xlm(find(w8)-1)).*xlm(w8)))*fs/(2*pi);
  fprintf('fd = %.1f MHz, fc = %.1f MHz: envelope error %.2f %%, residual rotation %.3f MHz (expect %.1f)\n', ...
          fdm, fcm, ephm, rot, f0 - fdm);
end

%% Figures: spectra at each stage, envelopes in time
F   = @(v) 20*log10(abs(fftshift(fft(v)))/max(abs(fft(x))) + 1e-12);
fax = ((0:N-1)' - N/2)*fs/N;
Hb  = 20*log10(abs(fftshift(fft(b, N))) + 1e-12);
figure('Position', [80 80 900 700]);
subplot(3,2,1); plot(fax, F(x));  axis([-20 20 -80 5]); title('1 RF (real): lobes at \pm f_0'); ylabel('dB');
subplot(3,2,2); plot(fax, F(xm)); axis([-20 20 -80 5]); title('2 mixed: baseband + image at -(f_0+f_d)');
subplot(3,2,3); plot(fax, F(xl/2), 'b', fax, Hb, 'k--'); axis([-20 20 -80 5]);
title(sprintf('3 after LPF (f_c = %g MHz, %d taps), filter dashed', fc, L)); ylabel('dB');
Nq = numel(iq); faq = ((0:Nq-1)' - Nq/2)*fsIQ/Nq;
subplot(3,2,4); plot(faq, 20*log10(abs(fftshift(fft(iq)))/max(abs(fft(iq))) + 1e-12));
axis([-fsIQ/2 fsIQ/2 -80 5]); title(sprintf('4 decimated by %d: span \\pm%.2f MHz', M, fsIQ/2)); xlabel('frequency (MHz)');
subplot(3,2,5); tw = t > 6 & t < 16;
plot(t(tw), abs(xa(tw)), 'k', t(idx(idx>=find(tw,1) & idx<=find(tw,1,'last'))), ...
     abs(iq(idx>=find(tw,1) & idx<=find(tw,1,'last'))), 'ro');
xlabel('time (\mus)'); ylabel('envelope'); title('|x_a| (k) and 2|IQ| at kept samples (o)');
subplot(3,2,6); plot(Ms, eM, 'bo-'); xlabel('decimation M'); ylabel('reconstruction error (%)');
title('band-limited reconstruction error vs M'); grid on;

figure('Position', [100 100 800 300]);
xlm = iqchain(x, fdm, 4);
plot(t(tw), real(xl(tw)), 'b', t(tw), real(xlm(tw)), 'r', t(tw), abs(xlm(tw)), 'k');
xlabel('time (\mus)'); ylabel('amplitude');
title('I with f_d = f_0 (b) and f_d = f_0 + 1 MHz (r): residual 1 MHz rotation, same envelope (k)');
