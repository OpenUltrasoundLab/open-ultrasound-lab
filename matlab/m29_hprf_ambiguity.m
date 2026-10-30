%% m29_hprf_ambiguity.m
% Purpose : High-PRF pulsed Doppler. When the PRF exceeds c/(2 z_gate), the
%           echo sampled at delay tau after the latest pulse also contains
%           echoes of earlier pulses from the depths z_m = (c/2)(tau + m T),
%           m = 0, 1, 2, ... (the extra sample volumes). The Nyquist velocity
%           grows with the PRF, but flow in any of those depths enters the
%           spectrum, and shallower ones are less attenuated.
%           Part 1: sample-volume depths, Nyquist velocity and transmit
%                   timing margin for a 12 cm gate at several PRFs.
%           Part 2: slow-time simulation of a fast target jet at the gate and
%                   a second vessel elsewhere: which one the spectrum shows.
% Change  : f0, zg (gate), PRF list, vessel depths and velocities,
%           attenuation, sample-volume length.
% Expect  : f0 = 2.5 MHz, gate 12 cm, theta = 0: conventional PRF limit
%           6.42 kHz, v_Nyq 0.99 m/s. PRF 10 kHz: tau = 55.8 us, sample
%           volumes at 4.3, 12.0 and 19.7 cm, v_Nyq 1.54 m/s; the 4.3 cm
%           volume is 19.3 dB less attenuated (0.5 dB/cm/MHz). A 1.3 m/s jet
%           aliases at 6 kHz and is shown correctly at 10 kHz; a vessel at
%           4.3 cm appears in the spectrum, one at 5.5 cm does not.
rng(3);

%% Parameters
c     = 1540;        % speed of sound (m/s)
f0    = 2.5e6;       % transmit frequency (Hz)
zg    = 0.12;        % gate depth set by the operator (m)
zdisp = 0.20;        % displayed (insonated) depth (m)
att   = 0.5;         % attenuation, dB/cm/MHz one way
Lsv   = 3e-3;        % sample-volume length (FWHM), m
tmin  = 5e-6;        % minimum margin between a transmit and the sampling instant (s)
prf0  = c/(2*zg);
fprintf('Gate %.0f cm, f0 %.1f MHz: conventional PRF <= %.2f kHz, v_Nyq = %.3f m/s (theta = 0)\n', zg*100, f0/1e6, prf0/1e3, c*prf0/(4*f0));

%% Part 1: sample-volume depths versus PRF
fprintf('Part 1:  PRF(kHz)  tau(us)  margin(us)  v_Nyq(m/s)  sample volumes within %.0f cm (cm)  [two-way attenuation vs gate, dB]\n', zdisp*100);
for PRF = [6 8 10 12 15 20]*1e3
    T = 1/PRF;  tau = mod(2*zg/c, T);  zm = c/2*(tau + (0:20)*T);  zm = zm(zm <= zdisp + 1e-9);
    marg = min(tau, T - tau);  rel = 2*att*f0/1e6*(zg - zm)*100;     % dB less attenuation than the gate
    s = sprintf('%5.1f ', zm*100);  r = sprintf('%+5.1f ', rel);
    flag = '';  if marg < tmin, flag = '  <- sampling too close to a transmit'; end
    fprintf('         %5.1f    %6.1f    %6.1f      %5.3f      %s [%s]%s\n', PRF/1e3, tau*1e6, marg*1e6, c*PRF/(4*f0), s, r, flag);
end

%% Part 2: what the spectrum shows
vj = 1.3;   vb = -0.4;          % jet at the gate (toward the probe), second vessel (away), m/s
Nfft = 256;  nseg = 60;  w = 0.5 - 0.5*cos(2*pi*(0:Nfft-1)'/Nfft);
cfg = [6e3 0.055; 10e3 0.043; 10e3 0.055];   % PRF, depth of the second vessel (m)
fprintf('Part 2: jet %.1f m/s at %.0f cm, second vessel %.1f m/s, equal backscatter, SV FWHM %.0f mm\n', vj, zg*100, vb, Lsv*1e3);
for ic = 1:size(cfg,1)
    PRF = cfg(ic,1);  zb = cfg(ic,2);  T = 1/PRF;  tau = mod(2*zg/c, T);  zm = c/2*(tau + (0:20)*T);
    gsv = @(z) sum(exp(-4*log(2)*(z - zm).^2/Lsv^2));              % sensitivity of the sum of sample volumes
    amp = [gsv(zg)*10^(-2*att*f0/1e6*zg*100/20), gsv(zb)*10^(-2*att*f0/1e6*zb*100/20)];
    fd  = 2*[vj vb]*f0/c;  sf = 0.03*abs(fd);                     % Doppler shifts and 3 % spectral width
    P = zeros(Nfft,1);  n = (0:Nfft-1)';
    for s = 1:nseg
        z = 0.002*amp(1)*(randn(Nfft,1) + 1j*randn(Nfft,1))/sqrt(2);   % noise 54 dB below the jet
        for k = 1:2
            ph = cumsum(2*pi*(fd(k) + sf(k)*randn(Nfft,1))/PRF);      % random-walk phase: narrowband process
            z = z + amp(k)*exp(1j*(ph + 2*pi*rand));
        end
        P = P + abs(fftshift(fft(z.*w))).^2;
    end
    fax = ((0:Nfft-1)' - Nfft/2)*PRF/Nfft;  vax = fax*c/(2*f0);
    [~, i1] = max(P);  Pd = 10*log10(P/max(P));
    vjw = mod(vj + c*PRF/(4*f0), c*PRF/(2*f0)) - c*PRF/(4*f0);  % where the jet lands (wrapped)
    [~, ij] = min(abs(vax - vjw));  [~, ib] = min(abs(vax - vb));
    fprintf('   PRF %4.1f kHz, vessel at %.1f cm: v_Nyq %.2f m/s; strongest line at %+.2f m/s; jet line %+.2f m/s at %5.1f dB; vessel line at %5.1f dB\n', ...
        PRF/1e3, zb*100, c*PRF/(4*f0), vax(i1), vjw, max(Pd(max(ij-2,1):min(ij+2,Nfft))), max(Pd(max(ib-2,1):min(ib+2,Nfft))));
end

%% Figure: sample-volume depths versus PRF
figure('Position', [50 50 800 380]); hold on;
prfs = linspace(4e3, 25e3, 400);
for m = 0:12
    zz = c/2*(mod(2*zg/c, 1./prfs) + m./prfs);  zz(zz > zdisp) = NaN;  plot(prfs/1e3, zz*100, 'b.', 'MarkerSize', 3);
end
plot([prf0 prf0]/1e3, [0 zdisp*100], 'k--'); plot([4 25], [zg zg]*100, 'r-');
xlabel('PRF (kHz)'); ylabel('sample-volume depth (cm)'); set(gca, 'YDir', 'reverse'); grid on;
title(sprintf('HPRF sample volumes for a %.0f cm gate (dashed: conventional limit)', zg*100));
