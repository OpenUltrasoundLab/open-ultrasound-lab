%% m04_rf_convolution.m  -- Module 04: RF line = pulse * reflectivity, spectra multiply, axial resolution
% Purpose : model an RF line as a Gaussian pulse convolved with a reflectivity profile,
%           verify that spectra multiply (with enough zero padding), measure the -6 dB
%           pulse length and axial resolution, resolve two reflectors, and compute the
%           envelope with an FFT-based analytic signal (no hilbert()).
% Change  : f0, fractional bandwidth FBW, the reflector spacings, fs.
% Expect  : tau(-6 dB) = 4 ln2/(pi B) = 0.294 us and dz = c tau/2 = 0.227 mm at 5 MHz, 60 %;
%           FFT and direct convolution agree to ~1e-13 when Nfft >= L+M-1;
%           Gabor product sigma_t' sigma_f' = 1/(4 pi) = 0.0796.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

%% Parameters
c   = 1.540;            % speed of sound (mm/us)
fs  = 100;              % sampling rate (MHz)
f0  = 5;                % center frequency (MHz)
FBW = 0.60;             % fractional bandwidth at -6 dB (amplitude)
B   = FBW*f0;                         % -6 dB bandwidth (MHz)
sig_t = sqrt(2*log(2))/(pi*B);        % envelope std (us): B = 2 sqrt(2 ln2) sigma_f, sigma_f = 1/(2 pi sigma_t)
tau_th = 4*log(2)/(pi*B);             % -6 dB pulse duration (us)
dz_th  = c*tau_th/2;                  % axial resolution (mm)

%% Pulse (two-way impulse response)
tp = (-ceil(4.5*sig_t*fs):ceil(4.5*sig_t*fs))/fs;     % time axis of the pulse (us)
p  = exp(-tp.^2/(2*sig_t^2)) .* cos(2*pi*f0*tp);
M  = numel(p);

%% Analytic signal by FFT (replaces hilbert): keep DC, double positive bins, zero negative bins
analytic = @(x, Nf) ifft( fft(x, Nf) .* [1, 2*ones(1, Nf/2-1), 1, zeros(1, Nf/2-1)] );

%% Measured pulse length from its envelope
Nf  = 2^nextpow2(8*M);
env_p = abs(analytic(p, Nf)); env_p = env_p(1:M);
above = find(env_p >= max(env_p)/2);
% sub-sample crossing by linear interpolation at both edges
i1 = above(1); i2 = above(end);
t1 = tp(i1-1) + (max(env_p)/2 - env_p(i1-1))/(env_p(i1) - env_p(i1-1))/fs;
t2 = tp(i2) + (env_p(i2) - max(env_p)/2)/(env_p(i2) - env_p(i2+1))/fs;
tau_meas = t2 - t1;
% Gabor product (rms widths of |envelope|^2 and |P(f)|^2, one-sided spectrum)
e2 = env_p.^2; mt = sum(tp.*e2)/sum(e2); st2 = sqrt(sum((tp-mt).^2.*e2)/sum(e2));
Pf = abs(fft(p, Nf)); fax = (0:Nf-1)*fs/Nf; pos = 1:Nf/2;
q2 = Pf(pos).^2; mf = sum(fax(pos).*q2)/sum(q2); sf2 = sqrt(sum((fax(pos)-mf).^2.*q2)/sum(q2));
fprintf('--- pulse: f0 = %.1f MHz, B(-6 dB) = %.2f MHz ---\n', f0, B);
fprintf('tau(-6 dB): theory 4ln2/(pi B) = %.4f us, measured from FFT envelope = %.4f us\n', tau_th, tau_meas);
fprintf('axial resolution c*tau/2: theory %.4f mm, measured %.4f mm\n', dz_th, c*tau_meas/2);
fprintf('Gabor product sigma_t''*sigma_f'' = %.5f (1/(4 pi) = %.5f)\n', st2*sf2, 1/(4*pi));

%% Reflectivity: speckle line (z = 0..6 mm) and the RF line by FFT convolution
zmax = 6; Ns = round(2*zmax/c*fs);              % samples covering 0..zmax (round trip)
rho = zeros(1, Ns); hit = rand(1, Ns) < 0.15;    % ~20 scatterers per mm (sample = 7.7 um)
rho(hit) = randn(1, nnz(hit));
L = Ns;
Nlin = L + M - 1;                                % linear convolution length
Nfft = 2^nextpow2(Nlin);
r_fft = real(ifft(fft(rho, Nfft) .* fft(p, Nfft)));  r_fft = r_fft(1:Nlin);
r_dir = conv(rho, p);
Nshort = L;                                      % too short (< L+M-1): circular wrap-around
r_circ = real(ifft(fft(rho, Nshort) .* fft(p, Nshort)));
fprintf('\n--- convolution ---\n');
fprintf('L = %d, M = %d, L+M-1 = %d, Nfft = %d: max|FFT conv - conv| = %.2e\n', L, M, Nlin, Nfft, max(abs(r_fft - r_dir)));
fprintf('Nfft = %d (< L+M-1): max error on first %d samples = %.2e (wrap-around)\n', Nshort, M, max(abs(r_circ(1:M) - r_dir(1:M))));
% spectra multiply
R = fft(r_dir, Nfft); S = fft(rho, Nfft); P = fft(p, Nfft);
fprintf('max|R - P.S| / max|R| = %.2e  (convolution in time = multiplication in frequency)\n', max(abs(R - P.*S))/max(abs(R)));

%% Two reflectors: axial resolution test
% The envelope of two echoes depends on their relative carrier phase phi = 4 pi f0 s / c:
% in phase (phi = 0) the envelopes add, in quadrature (phi = pi/2) the intensities add
% (= the average over random phase), in antiphase they cancel at the midpoint.
spac = [0.15 0.20 dz_th 0.25 0.30 0.40];           % spacings (mm)
fprintf('\n--- two equal reflectors at 3 mm: envelope dip between them (dB) ---\n');
fprintf('%10s %10s %12s %10s %10s %12s\n', 'spacing mm', 'spacing/dz', 'phi (deg)', 'dip RF', 'dip phi=0', 'dip phi=90');
zc = 3; Nr = 2^nextpow2(Ns + M);
gEnv = @(t) exp(-t.^2/(2*sig_t^2));               % Gaussian envelope
tt = linspace(-1, 1, 4001);                        % fine time axis (us) for the envelope-only cases
dipOf = @(e, a, b) 20*log10(min(e(a:b))/min(e(a), e(b)));
figure('Name', 'two reflectors', 'Position', [100 100 760 560]);
for is = 1:numel(spac)
    s = spac(is);
    rr = zeros(1, Ns);
    k1 = round(2*(zc - s/2)/c*fs) + 1; k2 = round(2*(zc + s/2)/c*fs) + 1;
    rr(k1) = 1; rr(k2) = rr(k2) + 1;
    rf = conv(rr, p);                                % pulse centered: shift by (M-1)/2
    rf = rf((M+1)/2 : (M+1)/2 + Ns - 1);
    en = abs(analytic(rf, Nr)); en = en(1:Ns);
    dip = dipOf(en, k1, k2);
    phi = mod(2*pi*f0*(k2 - k1)/fs, 2*pi);          % actual relative phase after sample rounding
    dt = s/c*2/2;                                    % half the echo separation (us)
    e1 = gEnv(tt + dt); e2 = gEnv(tt - dt);
    [~, a] = min(abs(tt + dt)); [~, b] = min(abs(tt - dt));
    d0 = dipOf(e1 + e2, a, b);                       % in phase
    d90 = dipOf(sqrt(e1.^2 + e2.^2), a, b);          % quadrature = intensities add
    fprintf('%10.3f %10.2f %12.0f %10.1f %10.1f %12.1f\n', s, s/dz_th, phi*180/pi, dip, d0, d90);
    if is == 2 || is == 5
        zz = (0:Ns-1)/fs*c/2;
        subplot(2, 1, 1 + (is == 5));
        plot(zz, rf, 'Color', [0.2 0.4 0.8]); hold on; plot(zz, en, 'r', 'LineWidth', 1.4);
        plot([zc-s/2 zc+s/2], [0 0], 'kv', 'MarkerFaceColor', 'k'); hold off;
        xlim([zc-1 zc+1]); xlabel('depth z = ct/2 (mm)'); ylabel('RF, envelope');
        title(sprintf('s = %.2f mm = %.2f dz, phase %.0f deg: dip %.1f dB', s, s/dz_th, phi*180/pi, dip));
    end
end
fprintf('At spacing = c tau/2 (%.3f mm) the quadrature (random-phase average) dip is -3.0 dB.\n', dz_th);

%% Figure: pulse, spectra, speckle RF line with FFT envelope
figure('Name', 'pulse and RF line', 'Position', [100 100 760 560]);
subplot(2, 2, 1);
plot(tp, p, 'b', tp, env_p, 'r', tp, -env_p, 'r'); grid on;
xlim([-0.6 0.6]); xlabel('time (us)'); ylabel('p(t)'); title(sprintf('pulse, tau = %.3f us', tau_meas));
subplot(2, 2, 2);
fa = (0:Nfft/2)*fs/Nfft;
Sn = abs(S(1:Nfft/2+1)); Pn = abs(P(1:Nfft/2+1)); Rn = abs(R(1:Nfft/2+1));
plot(fa, Sn/max(Sn), 'Color', [0.6 0.6 0.6]); hold on;
plot(fa, Pn/max(Pn), 'b', 'LineWidth', 1.5); plot(fa, Rn/max(Rn), 'r'); hold off;
xlim([0 12]); xlabel('frequency (MHz)'); ylabel('|.| (norm.)'); title('|S| gray, |P| blue, |PS| red');
subplot(2, 1, 2);
rfs = r_dir((M+1)/2 : (M+1)/2 + Ns - 1);
ens = abs(analytic(rfs, 2^nextpow2(2*Ns))); ens = ens(1:Ns);
zz = (0:Ns-1)/fs*c/2;
plot(zz, rfs, 'Color', [0.2 0.4 0.8]); hold on; plot(zz, ens, 'r', 'LineWidth', 1.2); hold off;
xlabel('depth z = ct/2 (mm)'); ylabel('amplitude'); title('speckle RF line and FFT-analytic envelope');
