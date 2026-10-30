%% m29_spectrogram_indices.m
% Purpose : Spectral (PW) Doppler from slow-time IQ to clinical numbers.
%           Part 1: STFT windows (rect, Hann, Blackman): -6 dB width,
%                   sidelobes, time/velocity resolution, and the window that
%                   balances resolution against the systolic chirp.
%           Part 2: Welch averaging inside a fixed window: std/mean of a
%                   noise-only bin versus the number of segments.
%           Part 3: synthetic pulsatile PW signal (parabolic profile, whole
%                   lumen, scatterers with finite lifetime, optional aperture
%                   broadening), Hann STFT, maximum-frequency envelope by a
%                   threshold and a percentile method, and the indices PSV,
%                   EDV, TAMAX, PI, RI, S/D, VTI against the true waveform.
%           Part 4: volume flow Q = TAMEAN x area and its error budget.
%           Part 5: simplified Bernoulli dP = 4 v^2 against 0.5 rho v^2.
% Change  : PRF, Nw, hop, SNR list, F-number, threshold, percentile, wall
%           filter cut-off, vessel diameter, angle error.
% Expect  : f0 = 5 MHz, theta = 60 deg, PRF = 10 kHz: v_Nyq = 1.54 m/s;
%           peak f_D = 3247 Hz. Hann -6 dB width 2.00 bins (rect 1.21,
%           Blackman 2.30); Nw = 128 -> 12.8 ms, 156 Hz, 0.048 m/s.
%           True waveform: PSV 1.000, EDV 0.203, TAMAX 0.276 m/s, PI 2.89
%           (EDV form), RI 0.797, S/D 4.94, VTI 27.6 cm. A threshold fixed
%           10 dB above the noise reads EDV and TAMAX high at 30 dB SNR (RI
%           low) and close to the truth at 10 dB SNR; at 3 dB it traces
%           noise. The 95th percentile reads PSV about 6 % low at any SNR.
%           Aperture broadening (F# 2) raises the threshold PSV by about 20 %
%           while RI barely moves; a wall filter above the diastolic Doppler
%           removes EDV (RI -> 1).
rng(2);

%% Parameters
c     = 1540;          % speed of sound (m/s)
f0    = 5e6;           % transmit frequency (Hz)
theta = 60*pi/180;     % beam-to-flow angle (rad)
PRF   = 10e3;          % pulse repetition frequency (Hz)
Tcyc  = 1;             % cardiac period (s), 60 beats per minute
kf    = 2*f0*cos(theta)/c;      % Doppler Hz per m/s along the vessel
fprintf('kf = %.1f Hz per m/s, v_Nyq = %.3f m/s, peak f_D = %.0f Hz\n', kf, c*PRF/(4*f0*cos(theta)), kf*1.0);

% periodic illustrative waveform of the centreline (maximum) velocity, m/s
gw   = @(t, m, s) exp(-((mod(t - m + 0.5, 1) - 0.5)/s).^2);
vfun = @(t) 0.2 + 0.8*gw(mod(t/Tcyc,1), 0.12, 0.05) - 0.08*gw(mod(t/Tcyc,1), 0.27, 0.025) + 0.1*gw(mod(t/Tcyc,1), 0.34, 0.05);
tt = (0:1e-5:1-1e-5)';  vv = vfun(tt);
tru.PSV = max(vv);  tru.Vmin = min(vv);  tru.EDV = vfun(0);  tru.TAMAX = mean(vv);
tru.PI = (tru.PSV - tru.EDV)/tru.TAMAX;  tru.PIg = (tru.PSV - tru.Vmin)/tru.TAMAX;
tru.RI = (tru.PSV - tru.EDV)/tru.PSV;  tru.SD = tru.PSV/tru.EDV;  tru.VTI = tru.TAMAX*Tcyc*100;
fprintf('True waveform: PSV %.3f, EDV %.3f, Vmin %.3f (notch), TAMAX %.3f m/s\n', tru.PSV, tru.EDV, tru.Vmin, tru.TAMAX);
fprintf('  PI (PSV-EDV)/TAMAX = %.2f, PI (max-min)/TAMAX = %.2f, RI = %.3f, S/D = %.2f, VTI = %.1f cm\n', ...
    tru.PI, tru.PIg, tru.RI, tru.SD, tru.VTI);

%% Part 1: windows and resolution
fprintf('Part 1: windows (N = 128, zero-padded x64)\n   window    -6dB width(bins)  peak sidelobe(dB)\n');
N = 128;  n = (0:N-1)';  Z = 64*N;
wins = {ones(N,1), 0.5 - 0.5*cos(2*pi*n/N), 0.42 - 0.5*cos(2*pi*n/N) + 0.08*cos(4*pi*n/N)};
wname = {'rect', 'Hann', 'Blackman'};
for k = 1:3
    W = abs(fft(wins{k}, Z));  W = W(1:Z/2)/W(1);
    i6 = find(W < 0.5, 1);  f6 = (i6 - 2 + (W(i6-1) - 0.5)/(W(i6-1) - W(i6)))/64;   % -6 dB half width (bins), interpolated
    knull = find(diff(W) > 0, 1);  w6(k) = 2*f6;
    fprintf('   %-9s %8.2f          %8.1f\n', wname{k}, 2*f6, 20*log10(max(W(knull:end))));
end
fprintf('   Hann at PRF %.0f kHz:  Nw   dt(ms)  df-6dB(Hz)  dv(m/s, 5 MHz, 60 deg)\n', PRF/1e3);
for Nw = [32 64 128 256]
    df = w6(2)*PRF/Nw;
    fprintf('                         %4d  %6.1f   %7.1f      %.3f\n', Nw, Nw/PRF*1e3, df, df/kf);
end
dvdt = max(diff(vv))/1e-5;  fdot = dvdt*kf;  Topt = sqrt(2/fdot);
fprintf('   max dv/dt = %.1f m/s^2 -> df/dt = %.1f kHz/s; smear |fdot|T + 2/T is smallest at T = %.1f ms (Nw = %.0f)\n', ...
    dvdt, fdot/1e3, Topt*1e3, Topt*PRF);

%% Part 2: Welch averaging inside one window (white noise only)
Nw = 128;  nwin = 3000;  fprintf('Part 2: Welch, Nw = %d, Hann segments with 50 %% overlap, %d noise windows\n   K   seg len   std/mean measured   predicted\n', Nw, nwin);
for K = [1 2 4 8]
    L = floor(2*Nw/(K+1));  hop = floor(L/2);  w = 0.5 - 0.5*cos(2*pi*(0:L-1)'/L);
    P = zeros(nwin, 1);
    for r = 1:nwin
        x = (randn(Nw,1) + 1j*randn(Nw,1))/sqrt(2);  s = 0;
        for q = 0:K-1
            X = fft(x(q*hop + (1:L)) .* w, 256);  s = s + abs(X(65))^2;   % one bin at PRF/4
        end
        P(r) = s/K;
    end
    rho = (sum(w(1:hop).*w(hop+1:2*hop))/sum(w.^2))^2;   % power correlation of 50 %-overlapped segments
    pred = sqrt((1 + 2*(1 - 1/K)*rho*(K > 1))/K);
    fprintf('   %d   %4d       %.3f               %.3f\n', K, L, std(P)/mean(P), pred);
end

%% Part 3: synthetic pulsatile PW signal, STFT, envelope, indices
t0 = -0.1;  Ns = round(2.2*PRF);  tn = t0 + (0:Ns-1)'/PRF;     % 2 cycles plus margins
fmax = kf*vfun(tn);  F = cumsum(fmax)/PRF;                      % cumulative phase per unit velocity fraction
Lf = 300;  Ksc = 64;  wl = 0.5 - 0.5*cos(2*pi*(0:Lf-1)'/Lf);    % 30 ms scatterer lifetime, 64 alive at once
NwS = 128;  hopS = 20;  wS = 0.5 - 0.5*cos(2*pi*(0:NwS-1)'/NwS);  nf = 256;
fax = ((0:nf-1)' - nf/2)*PRF/nf;  pos = find(fax > 0);  neg = find(fax < -300);
cols = 1:hopS:(Ns - NwS + 1);  tc = tn(cols + NwS/2);  ncol = numel(cols);
fprintf('Part 3: Hann Nw = %d (%.1f ms), hop %d (%.0f columns/s), %d scatterers, lifetime %.0f ms\n', NwS, NwS/PRF*1e3, hopS, PRF/hopS, Ksc, Lf/PRF*1e3);
cases = [Inf 30 10 0;  Inf 30 20 0;  Inf 10 10 0;  Inf 10 3 0;  2 30 10 0;  Inf 30 10 900];   % F#, SNR dB, thr dB, wall cut Hz
nc = size(cases, 1);  nreal = 3;  res = zeros(nc, 2, 6);  tamean = 0;   % indices averaged over realizations and cycles
mf = @(x) median([x([1 1:end-1]) x x([2:end end]) x([1 1 1:end-2]) x([3:end end end])], 2);   % 5-column median
for r = 1:nreal
    % same scatterers (velocity fraction, amplitude, beam angles) with and without aperture broadening
    zz = zeros(Ns, 2);  alpha2 = atan(1/(2*2));
    for k = 1:Ksc
        st = round((k-1)*Lf/Ksc) - Lf;
        while st < Ns
            idx = st + (1:Lf)';  ok = idx >= 1 & idx <= Ns;
            u = rand;  a = (randn + 1j*randn)/sqrt(2);  q1 = 2*rand - 1;  q2 = 2*rand - 1;
            b = [1, (cos(theta + alpha2*q1) + cos(theta + alpha2*q2))/(2*cos(theta))];
            dF = F(idx(ok)) - F(min(max(st+1,1),Ns));
            for ib = 1:2, zz(idx(ok),ib) = zz(idx(ok),ib) + a*wl(ok).*exp(2j*pi*u*b(ib)*dF); end
            st = st + Lf;
        end
    end
    wn = (randn(Ns,1) + 1j*randn(Ns,1))/sqrt(2);                    % unit-power noise shared by all cases
    for ic = 1:nc
        Fn = cases(ic,1);  snr = 10^(cases(ic,2)/10);  thr = 10^(cases(ic,3)/10);  fwall = cases(ic,4);
        z = zz(:, 1 + (Fn == 2));  z = z + sqrt(mean(abs(z).^2)/snr)*wn;
        S = zeros(nf, ncol);
        for j = 1:ncol, S(:,j) = abs(fftshift(fft(z(cols(j)-1+(1:NwS)).*wS, nf))).^2; end
        Pn = mean(mean(S(neg,:)));                                % noise floor from the reverse-flow half
        S(abs(fax) < fwall, :) = 0;                               % ideal wall filter
        envT = zeros(ncol,1);  envP = zeros(ncol,1);  fmean = zeros(ncol,1);
        for j = 1:ncol
            p = S(pos,j);  above = p > thr*Pn;
            run3 = above & [false; above(1:end-1)] & [false; false; above(1:end-2)];
            i3 = find(run3, 1, 'last');  if ~isempty(i3), envT(j) = fax(pos(i3)); end
            q = p - Pn;  q(fax(pos) < fwall) = 0;  cq = cumsum(q);
            if cq(end) > 0, envP(j) = fax(pos(find(cq >= 0.95*cq(end), 1))); end
            fmean(j) = max(sum(fax(pos).*q)/max(sum(q), eps), 0);
        end
        envs = {mf(envT)/kf, mf(envP)/kf};
        for im = 1:2
            e = envs{im};
            for cy = 0:1
                in = tc >= cy & tc < cy + 1;  [~, i0] = min(abs(tc - cy));
                PSV = max(e(in));  EDV = e(i0);  TAM = mean(e(in));
                res(ic,im,:) = res(ic,im,:) + reshape([PSV EDV TAM (PSV-EDV)/TAM (PSV-EDV)/PSV PSV/EDV], 1, 1, 6)/(2*nreal);
            end
        end
        if ic == 1
            tamean = tamean + mean(fmean(tc >= 0 & tc < 2))/kf/nreal;
            if r == 1, S1 = S;  e1 = envs; end
        end
    end
end
fprintf('   (%d realizations x 2 cycles; same scatterers and noise in every case)\n', nreal);
fprintf('   F#   SNR thr wall | method      PSV    EDV   TAMAX   PI    RI    S/D   VTI(cm)\n');
fprintf('   true              |           %.3f  %.3f  %.3f  %.2f  %.3f  %5.2f  %.1f\n', tru.PSV, tru.EDV, tru.TAMAX, tru.PI, tru.RI, tru.SD, tru.VTI);
mname = {'threshold', 'p95'};
for ic = 1:nc
    for im = 1:2
        R = squeeze(res(ic,im,:))';
        fprintf('   %-4g %3.0f %3.0f %4.0f | %-9s %6.3f  %.3f  %.3f  %.2f  %.3f  %5.2f  %.1f\n', cases(ic,1), cases(ic,2), cases(ic,3), cases(ic,4), mname{im}, R, R(3)*Tcyc*100);
    end
end
fprintf('   TAMEAN (intensity-weighted mean, first case) = %.3f m/s (TAMAX/2 = %.3f)\n', tamean, tru.TAMAX/2);

%% Part 4: volume flow and its error budget
d = 6e-3;  A = pi*d^2/4;  Q = tru.TAMAX/2*A*6e7;
fprintf('Part 4: d = %.1f mm, A = %.2f mm^2, TAMEAN = TAMAX/2 = %.3f m/s -> Q = %.0f mL/min\n', d*1e3, A*1e6, tru.TAMAX/2, Q);
dth = 5*pi/180;  dd = 0.3e-3;
ea = cos(theta)./cos(theta + [dth -dth]) - 1;  ed = ((d + [dd -dd])/d).^2 - 1;
fprintf('   angle 60 +/- 5 deg: %+.1f / %+.1f %% (linear %.1f %%); diameter %.1f +/- %.1f mm: %+.2f / %+.2f %% (linear %.1f %%); RSS %.1f %%\n', ...
    100*ea, 100*tan(theta)*dth, d*1e3, dd*1e3, 100*ed, 100*2*dd/d, 100*hypot(tan(theta)*dth, 2*dd/d));
fprintf('   profile: centreline-only sample of a parabolic profile reads TAMEAN = TAMAX -> Q x 2.00\n');
fprintf('   Womersley number (R = 3 mm, 1 Hz, rho 1060, mu 3.5 mPa s) = %.1f\n', d/2*sqrt(2*pi*1*1060/3.5e-3));

%% Part 5: simplified Bernoulli
rho = 1060;  mmHg = 133.322;
fprintf('Part 5: 0.5*rho = %.3f mmHg per (m/s)^2 (rule: 4)\n', 0.5*rho/mmHg);
for v = [2 3 4 5]
    fprintf('   v = %d m/s: 4v^2 = %3.0f mmHg, 0.5 rho v^2 = %5.1f mmHg; 20 deg misalignment -> %5.1f mmHg (%.1f %%)\n', ...
        v, 4*v^2, 0.5*rho*v^2/mmHg, 4*(v*cos(20*pi/180))^2, 100*(cos(20*pi/180)^2 - 1));
end

%% Figures
figure('Position', [50 50 1100 380]);
subplot(1,2,1); imagesc(tc, fax/kf, 10*log10(S1/max(S1(:)) + 1e-6)); axis xy; caxis([-50 0]); colormap(gray);
hold on; plot(tc, e1{1}, 'r', tc, e1{2}, 'c', tt, vv, 'y--'); xlabel('time (s)'); ylabel('velocity (m/s)');
title('Hann STFT, SNR 30 dB, envelopes'); legend('threshold 10 dB', 'p95', 'true v_{max}');
subplot(1,2,2); W = abs(fft(wins{2}, Z)); plot((0:Z/2-1)/64, 20*log10(W(1:Z/2)/W(1))); hold on;
W = abs(fft(wins{1}, Z)); plot((0:Z/2-1)/64, 20*log10(W(1:Z/2)/W(1))); xlim([0 8]); ylim([-80 0]); grid on;
xlabel('frequency (bins)'); ylabel('dB'); legend('Hann', 'rect'); title('Window spectra');
