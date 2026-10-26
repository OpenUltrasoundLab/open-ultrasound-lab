%% m08_attenuation_downshift.m  -  attenuation downshift and a toy QUS estimate
% Purpose : (1) multiply a Gaussian amplitude spectrum (center f0, std sigma_f)
%           by the round-trip attenuation exp(-2*alpha0*f*z) and compare the
%           measured power-spectrum centroid with f0 - 2*alpha0*sigma_f^2*z.
%           alpha0 is the AMPLITUDE attenuation slope in Np/cm/MHz, z is the
%           ONE-WAY depth in cm, sigma_f is the std of the AMPLITUDE spectrum
%           (= sqrt(2) x the std of the power spectrum).
%           (2) toy QUS: simulate speckle RF lines whose pulse downshifts with
%           depth, demodulate to IQ, estimate the local center frequency with
%           the fast-time lag-one autocorrelation, fit the slope versus depth
%           and invert it for alpha0.
% Change  : alpha_dB (dB/cm/MHz), f0, bw, window length, number of lines.
% Expect  : slope -0.187 MHz/cm at 5 MHz / 60 % / 0.5 dB/cm/MHz;
%           part (1) matches theory to < 0.001 MHz; part (2) recovers
%           alpha0 = 0.49 dB/cm/MHz (true 0.50; a few % spread from speckle
%           and from the LPF/decimation choices).
clear; close all; rng(1);

%% Parameters
f0 = 5;  bw = 0.6;                         % MHz, -6 dB fractional BW (amplitude)
sigf = bw*f0/(2*sqrt(2*log(2)));           % std of the amplitude spectrum (MHz)
alpha_dB = 0.5;                            % dB/cm/MHz (liver-like)
alpha = alpha_dB/(20*log10(exp(1)));       % Np/cm/MHz  (1 Np = 8.686 dB)
fprintf('sigma_f = %.3f MHz (power-spectrum sigma_p = %.3f MHz), alpha0 = %.4f Np/cm/MHz\n', ...
        sigf, sigf/sqrt(2), alpha);

%% 1. Spectral model: centroid versus depth
f  = linspace(0, 20, 4001)';               % MHz (positive frequencies only)
z  = 0:0.5:8;                              % one-way depth (cm)
G  = exp(-(f - f0).^2/(2*sigf^2));         % transmitted amplitude spectrum
S  = G .* exp(-2*alpha*f*z);               % received amplitude spectra (columns)
P  = S.^2;                                 % power spectra
fcP = sum(f.*P)./sum(P);                   % power-spectrum centroid (MHz)
fcA = sum(f.*S)./sum(S);                   % amplitude-spectrum centroid, for comparison
fcT = f0 - 2*alpha*sigf^2*z;               % theory
sdP = sqrt(sum((f - fcP).^2.*P)./sum(P));  % power-spectrum std (should stay sigma_f/sqrt 2)
fprintf('\nTheory slope -2*alpha0*sigma_f^2 = %.4f MHz/cm\n', -2*alpha*sigf^2);
fprintf(' z (cm)   fc power  fc ampl  theory   sigma_p\n');
for i = 1:2:numel(z)
  fprintf('%6.1f   %7.3f  %7.3f  %7.3f  %7.3f\n', z(i), fcP(i), fcA(i), fcT(i), sdP(i));
end
fprintf('max |fc power - theory| = %.2e MHz\n', max(abs(fcP - fcT)));
lossdB = 20*log10(exp(-2*alpha*f0*z(end)));
fprintf('Round-trip loss at f0, z = %.0f cm: %.1f dB (= 2 x %.1f dB/cm/MHz x %.0f MHz x %.0f cm)\n', ...
        z(end), lossdB, alpha_dB, f0, z(end));

%% 2. Toy QUS: speckle lines with a depth-dependent pulse
c = 1.54; fs = 40; M = 4; fsIQ = fs/M; fd = f0;      % mm/us, MHz
zmax = 80;  Nt = round(2*zmax/c*fs);  t = (0:Nt-1)'/fs;  zt = c*t/2;   % depth axis (mm)
Nl = 64;                                   % RF lines (independent speckle)
blk = 2;                                   % pulse is held constant over 2 mm blocks
sigt = 1/(2*pi*sigf);  tp = (-round(5*sigt*fs):round(5*sigt*fs))'/fs;
nb = ceil(zmax/blk);
RF = zeros(Nt, Nl);
for ib = 1:nb
  rows = find(zt >= (ib-1)*blk & zt < ib*blk);
  zc = ((ib-0.5)*blk)/10;                  % block center (cm)
  fcz = f0 - 2*alpha*sigf^2*zc;            % local center frequency (closed form)
  pulse = exp(-tp.^2/(2*sigt^2)).*cos(2*pi*fcz*tp);   % TGC: amplitude loss removed
  scat = zeros(Nt, Nl); scat(rows, :) = randn(numel(rows), Nl);
  RF = RF + conv2(scat, pulse, 'same');
end
% IQ chain: mix, 49-tap Hamming windowed-sinc LPF (fc = 4.5 MHz), decimate by M
L = 49; fcut = 4.5; k = (-(L-1)/2:(L-1)/2)'; u = 2*fcut/fs*k;
b = (2*fcut/fs)*(sin(pi*u) + (u==0))./(pi*u + (u==0)).*(0.54 - 0.46*cos(2*pi*(0:L-1)'/(L-1)));
b = b/sum(b);
IQ = 2*conv2(RF.*exp(-1j*2*pi*fd*t), b, 'same');  IQ = IQ(1:M:end, :);  zIQ = zt(1:M:end);
% Lag-one autocorrelation in 4 mm windows, averaged over lines
win = 4; zc = (win/2 + 5):2:(zmax - win/2 - 5);  fhat = zeros(size(zc));
for i = 1:numel(zc)
  r = find(abs(zIQ - zc(i)) < win/2);
  R1 = sum(sum(conj(IQ(r(1:end-1), :)).*IQ(r(2:end), :)));
  fhat(i) = fd + fsIQ/(2*pi)*angle(R1);
end
p = polyfit(zc/10, fhat, 1);               % slope in MHz/cm
alphaHat = -p(1)/(2*sigf^2);
fprintf('\nToy QUS (%d lines, %g mm windows, IQ at %.0f MHz):\n', Nl, win, fsIQ);
fprintf('  fitted slope = %.4f MHz/cm (theory %.4f)\n', p(1), -2*alpha*sigf^2);
fprintf('  alpha0 estimate = %.4f Np/cm/MHz = %.3f dB/cm/MHz (true %.3f)\n', ...
        alphaHat, alphaHat*20*log10(exp(1)), alpha_dB);
fprintf('  rms scatter of f_hat about the fit = %.3f MHz\n', std(fhat - polyval(p, zc/10)));

%% Figures
figure('Position', [80 80 950 360]);
subplot(1,3,1); plot(f, S(:, 1:4:end)./max(S(:,1))); xlim([0 12]);
xlabel('frequency (MHz)'); ylabel('amplitude (norm. to z = 0)');
title('received spectra, z = 0, 2, 4, 6, 8 cm');
subplot(1,3,2); plot(z, fcP, 'bo', z, fcT, 'k-', z, fcA, 'r.'); grid on;
xlabel('one-way depth z (cm)'); ylabel('center frequency (MHz)');
legend('power centroid', 'f_0 - 2\alpha_0\sigma_f^2 z', 'amplitude centroid', 'Location', 'southwest');
title('downshift: model check');
subplot(1,3,3); plot(zc/10, fhat, 'b.', zc/10, polyval(p, zc/10), 'r-', zc/10, f0 - 2*alpha*sigf^2*zc/10, 'k--');
grid on; xlabel('depth (cm)'); ylabel('f_c estimate (MHz)');
title(sprintf('toy QUS: \\alpha_0 = %.2f dB/cm/MHz', alphaHat*20*log10(exp(1))));
legend('R(1) estimate', 'linear fit', 'theory', 'Location', 'southwest');
