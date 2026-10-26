%% m02_leakage_zero_padding.m  -- Module 02: leakage, scalloping, Dirichlet kernel, zero padding
% Purpose : show that a tone between bins spreads over every bin (leakage), that the bins sample
%           the Dirichlet kernel D_N(delta) = sin(pi*delta)/sin(pi*delta/N), measure the scalloping
%           loss for rectangular and Hann windows, and show that zero padding samples the same
%           DTFT more finely without improving resolution.
% Change  : N, k0 (tone frequency in bins), the padding factor P, the tone separations sep.
% Expect  : complex tone, half-bin offset: rect -3.92 dB, Hann -1.42 dB; measured |X[k]|/N equals
%           |D_N(delta)|/N to ~1e-15; padded bins coincide with the DTFT; two tones 0.7 bin apart
%           show one peak at P = 1, 4, 16; two tones 2 bins apart show two peaks. Near 1/T the
%           outcome depends on the relative phase of the tones (section 4 sweeps it).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

N  = 64;                                   % record length (samples)
n  = (0:N-1).';
hann = 0.5 - 0.5*cos(2*pi*n/N);            % periodic Hann, written out (no toolbox)
dirichlet = @(d) (abs(d) < 1e-12) + (abs(d) >= 1e-12).*sin(pi*d)./(N*sin(pi*d/N) + (abs(d) < 1e-12));  % D_N(d)/N, = 1 at d = 0

%% 1. Scalloping loss: complex tone at an integer bin and half a bin off
fprintf('--- scalloping loss (complex tone, N = %d) ---\n', N);
for k0 = [10 10.25 10.5]
    z = exp(2j*pi*k0*n/N);
    Xr = fft(z);        Xh = fft(z.*hann);
    pr = max(abs(Xr))/N;                   % rect: on-bin peak = N
    ph = max(abs(Xh))/(N*mean(hann));      % Hann: on-bin peak = N * coherent gain (0.5)
    fprintf('k0 = %5.2f: peak bin rect %6.2f dB (theory %6.2f), Hann %6.2f dB\n', k0, 20*log10(pr), ...
            20*log10(abs(dirichlet(k0 - round(k0)))), 20*log10(ph));
end
fprintf('theory at half bin: 20log10(2/pi) = %.2f dB (rect)\n', 20*log10(2/pi));

%% 2. The bins sample the Dirichlet kernel
k0 = 10.3; z = exp(2j*pi*k0*n/N); X = fft(z);
k = (0:N-1).'; d = k0 - k;
fprintf('\n--- Dirichlet check, k0 = %.2f ---\n', k0);
fprintf('max | |X[k]|/N - |D_N(k0-k)|/N | = %.1e\n', max(abs(abs(X)/N - abs(dirichlet(d)))));
fprintf('bins 9..12: |X|/N = %s\n', num2str(abs(X(10:13)).'/N, '%8.4f'));
fprintf('fraction of energy outside bins 9..11: %.1f %%\n', 100*(1 - sum(abs(X(10:12)).^2)/sum(abs(X).^2)));

%% 3. Zero padding: the same DTFT on a finer grid
x = cos(2*pi*k0*n/N);                      % real tone this time
P = 8;
Xp = fft(x, P*N);                          % pads with (P-1)N zeros
fd = (0:P*N-1)/P;                          % frequency in bins of the original record
Xd = zeros(size(fd));                      % DTFT evaluated directly
for i = 1:numel(fd), Xd(i) = sum(x.*exp(-2j*pi*fd(i)*n/N)); end
fprintf('\n--- zero padding ---\nmax |padded FFT - direct DTFT| = %.1e\n', max(abs(Xp(:).' - Xd)));
[~, ip] = max(abs(Xp(1:P*N/2))); X1 = fft(x); [~, i1] = max(abs(X1(1:N/2)));
fprintf('padded peak at %.3f bins (true %.2f); unpadded peak at bin %d\n', fd(ip), k0, i1 - 1);

%% 4. Resolution is set by N, not by padding
fprintf('\n--- two tones: padding vs resolution (N = 16) ---\n');
Ns = 16; ns = (0:Ns-1).';
cnt = @(v) sum(v(2:end-1) > v(1:end-2) & v(2:end-1) > v(3:end) & v(2:end-1) > 0.5*max(v));   % strong local maxima
for sep = [0.5 2.0]
    xs = cos(2*pi*4*ns/Ns) + cos(2*pi*(4 + sep)*ns/Ns);
    msg = sprintf('separation %.1f bin:', sep);
    for Pp = [1 4 16]
        V = abs(fft(xs, Pp*Ns)); V = V(1:Pp*Ns/2 + 1);
        msg = [msg sprintf('  P = %2d -> %d peak(s)', Pp, cnt(V(:)))]; %#ok<AGROW>
    end
    fprintf('%s\n', msg);
end
fprintf('relative phase sweep (P = 16, 24 phases): fraction of phases with two peaks\n');
for sep = [0.5 1.0 1.5 2.0]
    nres = 0;
    for ph = (0:23)*2*pi/24
        xs = cos(2*pi*4*ns/Ns) + cos(2*pi*(4 + sep)*ns/Ns + ph);
        V = abs(fft(xs, 16*Ns)); V = V(1:8*Ns + 1); nres = nres + (cnt(V(:)) >= 2);
    end
    fprintf('  separation %.1f bin: %4.0f %%\n', sep, 100*nres/24);
end

%% 5. Figures
figure('Position', [100 100 900 600]);
subplot(2, 1, 1);
plot(fd, 20*log10(abs(Xd)/(N/2) + 1e-12), 'Color', [0.5 0.5 0.5]); hold on;
Xb = 20*log10(abs(X1(:))/(N/2) + 1e-12);
stem(0:N/2, Xb(1:N/2+1), 'filled', 'BaseValue', -80);
plot(fd(1:P*N/2), 20*log10(abs(Xp(1:P*N/2))/(N/2) + 1e-12), '.'); hold off;
xlim([0 N/2]); ylim([-80 5]); grid on;
xlabel('frequency (bins of the N-sample record)'); ylabel('dB re A N/2');
title(sprintf('Real tone at k_0 = %.1f: DTFT (grey), N-point bins (stems), %dx padded (dots)', k0, P));
subplot(2, 1, 2);
dd = linspace(-4, 4, 801);
plot(dd, 20*log10(abs(dirichlet(dd + 1e-9))), dd, 20*log10(abs(arrayfun(@(q) sum(hann.*exp(2j*pi*q*n/N)), dd))/(N/2)));
ylim([-80 3]); grid on; xlabel('offset \delta from the bin (bins)'); ylabel('dB');
legend('rectangular |D_N(\delta)|/N', 'Hann'); title('Window kernels: first sidelobe -13.3 dB vs -31.5 dB');
