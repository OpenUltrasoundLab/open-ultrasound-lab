%% m10_attenuation_mismatch.m
% A chirp through attenuating tissue, compressed with (a) the filter matched to
% the transmitted echo, (b) a depth-adapted matched filter, (c) an equalizing
% (Wiener-type) filter that restores the surface pulse shape.
%   - Attenuation alpha (dB/cm/MHz), linear in f, two-way path 2z, no dispersion.
%   - Received spectrum E(f) = S(f) H_tr(f) exp(-2 alpha_Np f z).
%   - SNR loss of (a) vs (b) = (int|E0||E|)^2 / (int|E0|^2 int|E|^2) <= 1 (Cauchy-Schwarz).
%   - With an exact (conjugate) filter the output spectrum is real and >= 0,
%     so the compressed envelope does not move; the raw echo's energy centroid
%     moves by about (fc - f0) T / B.
% What to change: alpha_dB, zlist, fbw, T, B, epsW.
% Expected output (alpha 0.5, 80 % transducer): centroid 5.00 -> 3.68 MHz at 6 cm,
% mismatch loss 2.0 dB at 6 cm and 4.6 dB at 10 cm; envelope shift 0.0 ns at
% every depth; at 6 cm the equalizer holds the -6 dB width at 0.38 us (0.34 at the
% surface, 0.48 for the adapted filter) but costs 5.8 dB.
% The adapted filter's PSL looks good only because its broad skirt has no null:
% read the -20 dB width (w20) as well. Runtime about 1 s.
clear; close all; rng(1);

%% Parameters
fs  = 50e6;                    % sampling rate (Hz)
f0  = 5e6;                     % chirp centre frequency (Hz)
T   = 10e-6;                   % chirp duration (s)
B   = 5e6;                     % swept bandwidth (Hz), 2.5-7.5 MHz
fbw = 0.8;                     % -6 dB fractional bandwidth of the two-way transducer response
alpha_dB = 0.5;                % attenuation (dB/cm/MHz)
zlist = [0 2 4 6 8 10];        % depth (cm)
epsW  = 1e-2;                  % equalizer regularization, relative to max|E|^2 (-20 dB)

Nfft = 4096; fax = [0:Nfft/2, -Nfft/2+1:-1]'*fs/Nfft;
n  = round(T*fs); t = ((0:n-1)' - (n-1)/2)/fs;
s  = cos(2*pi*(f0*t + B/(2*T)*t.^2));                 % RF up-chirp
S  = fft(s, Nfft);
sf = fbw*f0/2.3548;
Htr = exp(-(abs(fax)-f0).^2/(2*sf^2));                % Gaussian two-way transducer response
E0 = S.*Htr;                                          % echo spectrum at z = 0
an = zeros(Nfft,1); an(1) = 1; an(2:Nfft/2) = 2; an(Nfft/2+1) = 1;
env = @(Y) fftshift(abs(ifft(Y.*an)));                % envelope of a filter output, lag 0 centred
k0 = Nfft/2 + 1;
xR  = @(a,i,v) i-1 + (a(i-1)-v)./(a(i-1)-a(i));
xL  = @(a,i,v) i+1 - (a(i+1)-v)./(a(i+1)-a(i));
wf  = @(a,k,v) xR(a, k-1+find(a(k:end) < v, 1), v) - xL(a, k+1-find(a(k:-1:1) < v, 1), v);  % full width at level v
pslf = @(a,k) 20*log10(max(a([1:(k+1-find(diff(a(k:-1:1)) > 0, 1))-1, (k-1+find(diff(a(k:end)) > 0, 1))+1:numel(a)])));
pos = fax > 0; tt = (0:Nfft-1)'/fs;

fprintf('alpha = %.2f dB/cm/MHz, chirp %.1f-%.1f MHz, T = %.0f us, transducer %.0f%%\n', alpha_dB, (f0-B/2)/1e6, (f0+B/2)/1e6, T*1e6, fbw*100);
fprintf(' z(cm)  fc(MHz)  2*a*f0*z(dB) | fixed: w6 w20 (ns) PSL(dB) loss(dB) shift(ns) | adapted: w6 w20 PSL | equalizer: w6 w20 PSL loss | echo centroid shift (us) vs (fc-f0)T/B\n');
keep = cell(3,1);
for z = zlist
    E = E0.*exp(-alpha_dB/8.686*abs(fax)/1e6*2*z); % amplitude attenuation, alpha in Np/cm/MHz, path 2z
    P = abs(E).^2; fc = sum(fax(pos).*P(pos))/sum(P(pos));
    Hs = {E0, E, E.*abs(E0).^2./(abs(E).^2 + epsW*max(abs(E).^2))};
    w6 = zeros(1,3); w20 = w6; psl = w6; snr = w6; sh = w6;
    for k = 1:3
        e = env(E.*conj(Hs{k})); [pk, kp] = max(e); a = e/pk;
        w6(k) = wf(a,kp,0.5)/fs; w20(k) = wf(a,kp,0.1)/fs; psl(k) = pslf(a,kp); sh(k) = (kp - k0)/fs;
        snr(k) = pk^2/sum(abs(Hs{k}).^2);           % peak^2 / output noise power (common scale)
        if z == 6, keep{k} = a; end
    end
    ea = abs(ifft(E.*an)).^2; ea0 = abs(ifft(E0.*an)).^2;    % raw echo envelopes (no compression)
    dtc = sum(tt.*ea)/sum(ea) - sum(tt.*ea0)/sum(ea0);
    fprintf(' %4.0f   %5.2f     %5.1f      |  %4.0f %4.0f  %6.1f  %6.2f   %5.2f   | %4.0f %5.0f %6.1f | %4.0f %4.0f %6.1f %6.2f |  %6.2f vs %6.2f\n', ...
        z, fc/1e6, 2*alpha_dB*f0/1e6*z, w6(1)*1e9, w20(1)*1e9, psl(1), 10*log10(snr(1)/snr(2)), sh(1)*1e9, ...
        w6(2)*1e9, w20(2)*1e9, psl(2), w6(3)*1e9, w20(3)*1e9, psl(3), 10*log10(snr(3)/snr(2)), dtc*1e6, (fc-f0)*T/B*1e6);
    if z == 6, E6 = E; end
end
alpha0 = alpha_dB/8.686;       % Np/cm/MHz
fprintf('Gaussian-spectrum shift for reference: fc = f0 - 2 alpha0 sigma_f^2 z, sigma_f = %.2f MHz -> %.2f MHz at 6 cm (the flat chirp band makes the true shift smaller)\n', ...
    sf/1e6, f0/1e6 - 2*alpha0*(sf/1e6)^2*6);

%% Figures
figure('Position',[100 100 900 600]);
subplot(2,1,1);
plot(fax(pos)/1e6, abs(E0(pos))/max(abs(E0)), 'k', fax(pos)/1e6, abs(E6(pos))/max(abs(E6)), 'r', 'LineWidth', 1.2);
xlim([0 10]); grid on; xlabel('frequency (MHz)'); ylabel('|E(f)| (normalized)');
legend('z = 0', 'z = 6 cm', 'Location', 'northeast');
title(sprintf('Echo spectrum: chirp x transducer x attenuation (%.1f dB/cm/MHz)', alpha_dB));
subplot(2,1,2);
lag = ((1:Nfft)' - k0)/fs*1e6;
plot(lag, 20*log10(keep{1}+1e-9), 'k', lag, 20*log10(keep{2}+1e-9), 'b', lag, 20*log10(keep{3}+1e-9), 'r');
xlim([-12 12]); ylim([-70 3]); grid on; xlabel('lag (\mus)'); ylabel('envelope (dB)');
legend('fixed (matched to z = 0)', 'depth-adapted matched', 'equalizer to z = 0 shape', 'Location', 'northeast');
title('Compressed echo at 6 cm');
