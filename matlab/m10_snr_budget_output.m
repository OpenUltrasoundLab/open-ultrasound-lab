%% m10_snr_budget_output.m
% SNR budget of coded excitation and its acoustic-output cost.
%   - Penetration gain: Delta z = G / (2 alpha f0) with alpha in dB/cm/MHz,
%     two-way path, narrowband. Broadband check: SNR(z) = int |E0(f)|^2 exp(-4 alpha0 f z) df.
%   - Acoustic output at equal peak pressure: MI is unchanged; PII, ISPTA and
%     time-averaged power (hence TI) rise by the energy ratio, which is the
%     coding gain itself.
%   - At equal total energy, one coded firing and K averaged short pulses give
%     the same SNR (white noise): coding saves time and peak pressure, not energy.
% What to change: alpha_dB, f0, fbw, TB, PRF, p_r3 (derated peak pressure).
% Expected output: 15 dB at 0.5 dB/cm/MHz and 5 MHz -> 3.0 cm; a broadband
% (Gaussian, cut at -24 dB) spectrum gives 4.9 cm because it downshifts (2.8 MHz there);
% chirp TB 50 at MI 1.9 raises ISPTA.3 at a point (200 pulses/s) from 27 to
% 2003 mW/cm^2 (x75.2 = 18.8 dB, the coding gain; limit 720); one coded firing
% vs 75 averaged short pulses: SNR equal within Monte Carlo error (about 0.2 dB).
% Runtime about 3 s.
clear; close all; rng(1);

c = 1540;            % sound speed (m/s)
rhoc = 1.5e6;        % characteristic impedance of tissue (Rayl), for intensity p^2/(rho c)

%% 1. Penetration gain, narrowband rule
codes = {'chirp TB 50, rect', 'chirp TB 50, Hann on receive', 'Barker-13', 'Golay 16 (pair, 2 firings)'};
G = [10*log10(1.505*50), 10*log10(1.505*50) - 1.76, 10*log10(13), 10*log10(2*16)];
fprintf('Penetration gain Delta z = G/(2 alpha f0)  (cm)\n%-30s  G (dB)   a=0.5,f0=3   a=0.5,f0=5   a=0.7,f0=5   a=0.5,f0=7.5\n', 'code');
for k = 1:numel(codes)
    fprintf('%-30s  %5.2f     %6.2f       %6.2f       %6.2f       %6.2f\n', codes{k}, G(k), ...
        G(k)/(2*0.5*3), G(k)/(2*0.5*5), G(k)/(2*0.7*5), G(k)/(2*0.5*7.5));
end
fprintf('Reference: 15 dB at 0.5 dB/cm/MHz, 5 MHz -> %.2f cm\n', 15/(2*0.5*5));

%% 2. Broadband check: the attenuation slope flattens with depth because the spectrum downshifts
alpha_dB = 0.5; f0 = 5; fbw = 0.6;                    % dB/cm/MHz, MHz, -6 dB fractional bandwidth
sf = fbw*f0/2.3548; f = linspace(0.01, 12, 1200)';   % MHz
E0sq = exp(-(f-f0).^2/sf^2).*(abs(f-f0) <= fbw*f0);  % |E0(f)|^2: Gaussian amplitude spectrum cut at its -24 dB points
a0 = alpha_dB/8.686;                                  % Np/cm/MHz
z = linspace(0, 30, 601);
SNRdB = 10*log10(sum(repmat(E0sq,1,numel(z)).*exp(-4*a0*f*z), 1)/sum(E0sq));   % relative to z = 0
SNRnb = -2*alpha_dB*f0*z;                             % narrowband rule
budget = 40;                                          % short pulse: detectable down to -40 dB re z = 0
zsp_nb = budget/(2*alpha_dB*f0); zsp_bb = interp1(-SNRdB, z, budget);
fprintf('\nBroadband check (alpha %.1f, f0 %.0f MHz, %.0f%% band), short-pulse budget %d dB\n', alpha_dB, f0, fbw*100, budget);
fprintf('   G (dB)   depth gain narrowband   depth gain broadband   centre frequency there (MHz)\n');
for Gk = [10 15 20]
    zG = interp1(-SNRdB, z, budget+Gk);
    fprintf('   %4.0f        %5.2f cm              %5.2f cm               %4.2f\n', Gk, Gk/(2*alpha_dB*f0), zG - zsp_bb, sum(f.*E0sq.*exp(-4*a0*f*zG))/sum(E0sq.*exp(-4*a0*f*zG)));
end
fprintf('   short-pulse depth: narrowband %.2f cm, broadband %.2f cm\n', zsp_nb, zsp_bb);

%% 3. Acoustic output at equal derated peak pressure
fc  = 5e6; B = fbw*fc;                                % Hz
p_r3 = 1.9*sqrt(fc/1e6)*1e6;                          % derated peak rarefactional pressure at MI = 1.9 (Pa)
fs  = 200e6;
PRF = 200;           % effective pulse rate at one point (Hz): scanned B-mode, e.g. 40 frames/s x 5 overlapping lines
sft = B/2.3548; st = 1/(2*pi*sft);
tg  = (-6*st:1/fs:6*st)'; pg = -p_r3*exp(-tg.^2/(2*st^2)).*cos(2*pi*fc*tg);     % short pulse, rarefactional peak p_r3
T   = 50/B; tc = (-T/2:1/fs:T/2)'; pc = p_r3*cos(2*pi*(fc*tc + B/(2*T)*tc.^2)); % chirp TB 50
PII = @(p) sum(p.^2)/fs/rhoc;                         % pulse-intensity integral (J/m^2)
MI  = @(p) max(-p)/1e6/sqrt(fc/1e6);
Ig = PRF*PII(pg)/10; Ic = PRF*PII(pc)/10;             % W/m^2 -> mW/cm^2 is x0.1
fprintf('\nAcoustic output at one point, effective pulse rate %.0f Hz, derated peak %.2f MPa (MI 1.9 at %.0f MHz)\n', PRF, p_r3/1e6, fc/1e6);
fprintf('   short pulse: MI %.2f, ISPTA.3 %.0f mW/cm^2\n', MI(pg), Ig);
fprintf('   chirp TB 50: MI %.2f, ISPTA.3 %.0f mW/cm^2  (ratio %.1f = %.2f dB; coding gain %.2f dB)\n', MI(pc), Ic, Ic/Ig, 10*log10(Ic/Ig), 10*log10(1.505*50));
fprintf('   to meet ISPTA.3 <= 720 mW/cm^2: pulse rate <= %.0f Hz, or peak pressure x %.2f (MI %.2f)\n', PRF*720/Ic, min(1,sqrt(720/Ic)), MI(pc)*min(1,sqrt(720/Ic)));
fprintf('   same pulses at 5 kHz on every point (unfocused ultrafast-like, same peak): short %.0f, chirp %.0f mW/cm^2\n', Ig*5e3/PRF, Ic*5e3/PRF);
fprintf('   TI scales with time-averaged power: x %.1f at the same PRF\n', Ic/Ig);
fprintf('   max PRF set by depth: 5 cm %.0f Hz, 10 cm %.0f Hz, 15 cm %.0f Hz\n', c/(2*0.05), c/(2*0.10), c/(2*0.15));
for Tk = [5 10 20 40]
    fprintf('   blind zone c T / 2 for T = %2d us: %.1f mm\n', Tk, c*Tk*1e-6/2*1e3);
end

%% 4. Equal energy, equal SNR: one coded firing vs K averaged short pulses
K = round(1.505*50);                                  % energy ratio of the TB-50 chirp to the short pulse
Ntrial = 2000; sigma = 4;                             % complex white noise
fs = 25e6;                                            % complex baseband rate is enough here (B = 3 MHz)
tb = (-6*st:1/fs:6*st)'; g = exp(-tb.^2/(2*st^2));    % baseband short pulse, peak 1
tcb = (-T/2:1/fs:T/2)'; sc = exp(1j*pi*(B/T)*tcb.^2); % baseband chirp, peak 1
SNR = @(v) abs(mean(v))^2/mean(abs(v-mean(v)).^2);
yc = sc'*(repmat(sc,1,Ntrial) + sigma*(randn(numel(sc),Ntrial)+1j*randn(numel(sc),Ntrial))/sqrt(2));
yg = zeros(1,Ntrial);
for k = 1:K
    yg = yg + g'*(repmat(g,1,Ntrial) + sigma*(randn(numel(g),Ntrial)+1j*randn(numel(g),Ntrial))/sqrt(2));
end
fprintf('\nEqual energy: coded single firing %.2f dB, %d averaged short pulses %.2f dB (energy ratio %.2f dB)\n', ...
    10*log10(SNR(yc)), K, 10*log10(SNR(yg)), 10*log10(sum(abs(sc).^2)/sum(g.^2)));

%% Figures
figure('Position',[100 100 900 420]);
plot(z, SNRnb, 'k--', z, SNRdB, 'b', z, SNRdB + 15, 'r', 'LineWidth', 1.2); hold on;
plot([0 20], -budget*[1 1], 'Color', [0.5 0.5 0.5]);
ylim([-120 20]); grid on; xlabel('depth (cm)'); ylabel('echo SNR re short pulse at z = 0 (dB)');
legend('narrowband -2\alpha f_0 z', 'broadband short pulse', 'coded, +15 dB', 'detection threshold', 'Location', 'southwest');
title(sprintf('SNR budget: %.1f dB/cm/MHz, %.0f MHz, %.0f%% band', alpha_dB, f0, fbw*100));
