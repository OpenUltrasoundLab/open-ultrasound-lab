%% m29_doppler_spectrum_profile.m
% Purpose : Doppler power spectrum of flow in a straight vessel. The spectrum
%           is the histogram of f_D = 2 v f0 cos(theta)/c over the scatterers
%           in the sample volume, broadened by the spread of beam angles
%           across the aperture (geometric = transit-time broadening).
%           Part 1: plug, parabolic and blunted profiles, full-lumen gate.
%           Part 2: gate size and elevation beam (thin vs wide).
%           Part 3: geometric broadening and the max-velocity overestimate.
%           Part 4: simulated slow-time signal + averaged periodogram check.
%           Part 5: STFT window length vs frequency/time resolution.
% Change  : theta, vmax, a, F-number (alpha), gate length, Nw.
% Expect  : f0 = 5 MHz, theta = 60 deg, vmax = 0.8 m/s: fDmax = 2597 Hz.
%           Parabolic profile, whole lumen, wide elevation beam: flat
%           spectrum, mean = fDmax/2 = 1299 Hz. Thin elevation beam: the
%           spectrum rises toward fDmax. F# 2 (alpha = 14.0 deg): the
%           broadened spectrum extends to fDmax*(1 + alpha*tan(theta)) =
%           1.42 fDmax; at -20 dB of its plateau it ends near 1.33 fDmax.
%           Hann window, PRF 5 kHz: Nw = 128 -> bin 39.1 Hz, -6 dB width
%           2.00 bins = 78 Hz, 25.6 ms.
rng(1);

%% Parameters
c     = 1540;          % speed of sound (m/s)
f0    = 5e6;           % centre frequency (Hz)
theta = 60*pi/180;     % beam-to-flow angle (rad)
a     = 3e-3;          % vessel radius (m)
vmax  = 0.8;           % centreline velocity (m/s)
Fnum  = 2;             % F-number of tx and rx apertures
alpha = atan(1/(2*Fnum));   % aperture half-angle (rad)
fDmax = 2*vmax*f0*cos(theta)/c;
fprintf('fDmax = %.0f Hz (v = %.2f m/s, f0 = %.1f MHz, theta = %d deg)\n', fDmax, vmax, f0/1e6, round(theta*180/pi));

nb   = 400;                          % frequency bins over [-0.1, 1.6] fDmax
fax  = linspace(-0.1, 1.6, nb) * fDmax;  df = fax(2) - fax(1);
hist1 = @(f, w) accumarray(min(max(round((f(:) - fax(1))/df) + 1, 1), nb), w(:), [nb 1])' / (sum(w(:)) * df);

% Cross-section coordinates: u = in-plane distance from the axis (u = s sin(theta),
% s along the beam), y = elevation. Uniform fine grid over the square, masked to the lumen.
ng = 801;  [U, Y] = meshgrid(linspace(-a, a, ng));  R2 = (U.^2 + Y.^2) / a^2;  in = R2 < 1;
prof = {@(r2) ones(size(r2)), @(r2) 1 - r2, @(r2) 1 - r2.^2};   % plug, parabolic, blunted (n = 4)
pname = {'plug', 'parabolic', 'blunted n=4'};

%% Part 1: profile shape, full-lumen gate, wide elevation beam
S1 = zeros(3, nb);
fprintf('Part 1 (whole lumen, wide elevation beam):\n');
for k = 1:3
    v = vmax * prof{k}(R2(in));  f = 2*v*f0*cos(theta)/c;
    S1(k, :) = hist1(f, ones(size(f)));
    fprintf('  %-12s mean f = %6.0f Hz (%.3f fDmax)   mean v = %.3f m/s\n', pname{k}, mean(f), mean(f)/fDmax, mean(v));
end

%% Part 2: gate length and elevation beam (parabolic profile)
fprintf('Part 2 (parabolic):\n');
gates = [1 0.3];  S2 = zeros(4, nb);  ii = 0;
for elev = [1 0]                      % 1 = wide elevation beam, 0 = thin (y = 0 only)
    for g = gates                     % gate covers |u| < g*a (centred)
        ii = ii + 1;
        if elev, m = in & abs(U) < g*a; uu = U(m); yy = Y(m);
        else,    uu = linspace(-g*a, g*a, 4001)'; yy = zeros(size(uu)); end
        f = 2*vmax*(1 - (uu.^2 + yy.^2)/a^2)*f0*cos(theta)/c;
        S2(ii, :) = hist1(f, ones(size(f)));
        enames = {'thin', 'wide'};
        fprintf('  elevation %-5s gate %3.0f %% of lumen: mean f = %.3f fDmax, min f = %.3f fDmax\n', ...
            enames{elev + 1}, 100*g, mean(f)/fDmax, min(f)/fDmax);
    end
end

%% Part 3: geometric broadening (tx and rx angles uniform in theta +/- alpha)
% f = v f0 (cos(theta_t) + cos(theta_r))/c ; linearized: triangle of half-width alpha*tan(theta)*f_D
v = vmax * prof{2}(R2(in));  v = v(randperm(numel(v), 20000));   % random subset of the lumen
fD = 2*v*f0*cos(theta)/c;
nd = 41;  d = linspace(-alpha, alpha, nd);  [Dt, Dr] = meshgrid(d);
fac = (cos(theta + Dt(:)) + cos(theta + Dr(:))) / (2*cos(theta));     % exact angle factor
fB = fD(:) * fac';                                          % every scatterer x angle pair
S3 = hist1(fB(:), ones(numel(fB), 1));
plateau = median(S1(2, fax > 0.2*fDmax & fax < 0.8*fDmax));
iEnd = find(S3 > 0.01*plateau, 1, 'last');
fprintf('Part 3: alpha = %.1f deg (F# %d): full extent = %.3f fDmax (exact cos(theta-alpha)/cos(theta) = %.3f, linear 1 + alpha tan(theta) = %.3f)\n', ...
    alpha*180/pi, Fnum, max(fB(:))/fDmax, cos(theta - alpha)/cos(theta), 1 + alpha*tan(theta));
fprintf('        -20 dB envelope edge = %.3f fDmax -> apparent vmax = %.3f m/s (true %.2f)\n', fax(iEnd)/fDmax, vmax*fax(iEnd)/fDmax, vmax);
fprintf('        mean f unchanged? %.3f vs %.3f fDmax\n', mean(fB(:))/fDmax, mean(fD)/fDmax);

%% Part 4: simulated slow-time signal and averaged periodogram (parabolic, broadened)
% Each segment is an independent realization of K scatterers (uniform over the disk,
% random tx/rx angles within +/- alpha, complex Gaussian amplitudes) plus white noise.
PRF = 8e3;  Nw = 256;  nseg = 200;  K = 300;  t = (0:Nw-1)/PRF;
win = 0.5 - 0.5*cos(2*pi*(0:Nw-1)'/Nw);  P = zeros(Nw, 1);  fsum = 0;
for s = 1:nseg
    rr = sqrt(rand(K,1));
    fk = 2*vmax*(1 - rr.^2)*f0*cos(theta)/c .* ((cos(theta + alpha*(2*rand(K,1)-1)) + cos(theta + alpha*(2*rand(K,1)-1))) / (2*cos(theta)));
    amp = (randn(1,K) + 1j*randn(1,K))/sqrt(2);
    seg = (amp * exp(2j*pi*fk*t)).';                         % Nw x 1 slow-time IQ
    seg = seg + sqrt(K*1e-3) * (randn(Nw,1) + 1j*randn(Nw,1))/sqrt(2);   % SNR 30 dB
    P = P + abs(fft(seg .* win)).^2;  fsum = fsum + mean(fk);
end
P = fftshift(P);  fP = (-Nw/2:Nw/2-1)'*PRF/Nw;
fprintf('Part 4: periodogram mean f = %.0f Hz; model mean = %.0f Hz (K = %d scatterers x %d segments)\n', ...
    sum(fP.*P)/sum(P), mean(fB(:)), K, nseg);

%% Part 5: STFT resolution (Hann window)
fprintf('Part 5: STFT at PRF = 5 kHz, Hann window\n   Nw   bin(Hz)  -6dB width(Hz)  dt(ms)  dv at theta=0 (mm/s)\n');
PRF5 = 5e3;
for Nw = [32 64 128 256]
    w = 0.5 - 0.5*cos(2*pi*(0:Nw-1)'/Nw);  Z = 64*Nw;
    W = abs(fft(w, Z));  W = W / W(1);  k6 = find(W < 0.5, 1) - 1;   % first bin below -6 dB
    w6 = 2 * k6 / Z * PRF5;                                            % two-sided -6 dB width (Hz)
    fprintf('  %4d   %6.1f     %6.1f        %5.1f     %6.1f\n', Nw, PRF5/Nw, w6, Nw/PRF5*1e3, c*w6/(2*f0)*1e3);
end

%% Figures
figure('Position', [50 50 1100 380]);
subplot(1,3,1); plot(fax/fDmax, S1*fDmax, 'LineWidth', 1.5); xlabel('f / f_{D,max}'); ylabel('normalized power density');
title('Profile shape (whole lumen)'); legend(pname, 'Location', 'north'); ylim([0 5]); xlim([-0.1 1.6]); grid on;
subplot(1,3,2); plot(fax/fDmax, S2*fDmax, 'LineWidth', 1.5); xlabel('f / f_{D,max}'); ylim([0 8]); xlim([-0.1 1.6]); grid on;
title('Parabolic: gate and elevation beam'); legend('wide, 100%', 'wide, 30%', 'thin, 100%', 'thin, 30%', 'Location', 'northwest');
subplot(1,3,3); plot(fax/fDmax, S1(2,:)*fDmax, fax/fDmax, S3*fDmax, 'LineWidth', 1.5); hold on;
pl = median(P(fP > 0.2*fDmax & fP < 0.6*fDmax));  pm = median(S3(fax > 0.2*fDmax & fax < 0.6*fDmax));
plot(fP/fDmax, P/pl*pm*fDmax, 'k.-'); xlabel('f / f_{D,max}'); grid on; xlim([-0.1 1.6]);
title(sprintf('Broadening, F#%d, \\theta = %d deg', Fnum, round(theta*180/pi))); legend('no broadening', 'model', 'simulated periodogram');
