%% m32_global_sos_focus.m
% Global sound speed by focus-quality maximization in a two-layer medium.
% Purpose : Simulate channel data (0 deg plane-wave transmit, 64-element
%           receive, 5 MHz) from wire targets under a 10 mm fat layer
%           (1450 m/s) over liver (1580 m/s). Receive times are ray-traced
%           with Snell's law (refraction included). For each target, sweep
%           the assumed speed c_hat, beamform a small ROI tied to two-way
%           time, and keep the c_hat that maximizes the coherence factor.
%           Compare with two averages of the true profile:
%             harmonic (time-average) V_h = z / T,   T = sum d_i/c_i
%             RMS (Dix)                V_rms = sqrt( sum c_i^2 dT_i / T )
%           Then invert the measured averages to interval (local) speeds
%           with the straight-ray rule (dz/dT) and with Dix (d(V^2 T)/dT).
% Change  : c_fat, d_fat, c_liv, f0, the receive aperture Nel.
% Expect  : best c_hat follows V_rms within about 1 m/s (0.6-1.0 m/s
%           high) and sits 1.5-2.3 m/s above V_h; V_rms exceeds V_h by
%           w1 w2 (c2 - c1)^2 / (2 V_h), about 1.1 m/s at 40 mm; the interval
%           speed below the fat is recovered to < 1 m/s by Dix and to
%           0.5-2.5 m/s by the straight-ray rule; the shallow fat target
%           gives 1450.
% Runtime : about 3 s in Octave.
clear; close all; rng(1);

c_fat = 1450; d_fat = 10e-3; c_liv = 1580;   % layer speeds (m/s), fat thickness (m)
f0 = 5e6; fs = 40e6; sig_t = 0.12e-6;        % centre frequency, sampling, Gaussian envelope std (s)
Nel = 64; pitch = 0.3e-3;                    % receive elements, pitch (m)
xe = ((0:Nel-1) - (Nel-1)/2)*pitch;          % element positions (m)
zt = [6 15 20 25 30 35 40 45 50]*1e-3;       % wire depths on axis (m)
nT = numel(zt);

% ---- ray tracing: one-way time from element (xe,0) to wire (0,z), Snell refraction
Trec = zeros(nT, Nel); Ttx = zeros(nT, 1);
for k = 1:nT
    z = zt(k);
    if z <= d_fat
        Trec(k, :) = sqrt(z^2 + xe.^2)/c_fat; Ttx(k) = z/c_fat; continue;
    end
    d2 = z - d_fat; X = abs(xe);
    lo = zeros(size(X)); hi = (1/c_liv)*(1 - 1e-12)*ones(size(X));   % ray parameter p (s/m)
    for it = 1:60
        p = (lo + hi)/2;
        f = d_fat*p*c_fat./sqrt(1 - (p*c_fat).^2) + d2*p*c_liv./sqrt(1 - (p*c_liv).^2);
        hi(f > X) = p(f > X); lo(f <= X) = p(f <= X);
    end
    p = (lo + hi)/2;
    Trec(k, :) = d_fat./(c_fat*sqrt(1 - (p*c_fat).^2)) + d2./(c_liv*sqrt(1 - (p*c_liv).^2));
    Ttx(k) = d_fat/c_fat + d2/c_liv;         % 0 deg plane wave: vertical path through flat layers
end

% ---- channel data as baseband IQ (analytic RF times exp(-j w0 t))
t = (0:round(80e-6*fs))/fs;                  % fast time (s)
IQ = zeros(numel(t), Nel);
for k = 1:nT
    tau = Ttx(k) + Trec(k, :);               % two-way arrival per element (s)
    D = bsxfun(@minus, t(:), tau);           % time relative to arrival
    IQ = IQ + exp(-D.^2/(2*sig_t^2)).*exp(-1j*2*pi*f0*repmat(tau, numel(t), 1));
end
IQ = IQ + 0.01*(randn(size(IQ)) + 1j*randn(size(IQ)))/sqrt(2);   % -40 dB noise

% ---- sweep c_hat, coherence factor over a time-tied ROI around each wire
chat = 1400:1:1650;
CF = zeros(nT, numel(chat));
Nt = numel(t);
xr = (-0.1:0.05:0.1)*1e-3;                   % ROI lateral positions (m), inside the PSF main lobe
for k = 1:nT
    t0 = Ttx(k) + min(Trec(k, :));           % two-way apex time of the wire (s)
    tr = t0 + (-0.05:0.01:0.05)*1e-6;         % ROI in two-way time (s), +-39 um
    [TR, XR] = ndgrid(tr, xr); TR = TR(:); XR = XR(:);   % pixels, column vectors
    for q = 1:numel(chat)
        ch = chat(q); zh = ch*TR/2;          % pixel depth tied to time
        tau = zh/ch + sqrt(zh.^2 + bsxfun(@minus, XR, xe).^2)/ch;   % 0 deg tx + rx delay, Npix x Nel
        u = tau*fs + 1; i0 = floor(u); fr = u - i0; i0 = min(max(i0, 1), Nt - 1);
        lin = bsxfun(@plus, i0, (0:Nel-1)*Nt);                    % linear indices into IQ
        s = (1 - fr).*IQ(lin) + fr.*IQ(lin + 1);                    % linear interpolation of baseband IQ
        s = s.*exp(1j*2*pi*f0*tau);          % remodulate to the focused phase
        CF(k, q) = sum(abs(sum(s, 2)).^2)/(Nel*sum(abs(s(:)).^2));  % energy-weighted CF over the ROI
    end
end

% ---- results versus the two averages
Vh = zeros(nT, 1); Vr = zeros(nT, 1); best = zeros(nT, 1); T1 = zeros(nT, 1);
fprintf('  z (mm)   best c_hat   V_h (harm.)   V_rms (Dix)   best - V_h   best - V_rms   CF max\n');
for k = 1:nT
    z = zt(k); d1 = min(z, d_fat); d2 = max(0, z - d_fat);
    dT = [d1/c_fat, d2/c_liv]; T1(k) = sum(dT);
    Vh(k) = z/T1(k); Vr(k) = sqrt((c_fat^2*dT(1) + c_liv^2*dT(2))/T1(k));
    [m, i] = max(CF(k, :)); i = min(max(i, 2), numel(chat) - 1);
    y = CF(k, i-1:i+1); best(k) = chat(i) + 0.5*(y(1) - y(3))/(y(1) - 2*y(2) + y(3));   % parabolic refinement
    fprintf('  %5.0f   %9.1f   %10.1f   %10.1f   %9.1f   %11.1f   %6.3f\n', z*1e3, best(k), Vh(k), Vr(k), best(k) - Vh(k), best(k) - Vr(k), m);
end

% ---- interval speeds from the measured averages (one-way times T1 known from the echo times)
fprintf('\nInterval speed between successive wires (true liver %d m/s):\n', c_liv);
fprintf('  interval (mm)   straight-ray dz/dT   Dix sqrt(d(V^2 T)/dT)\n');
for k = 3:nT
    zs_k = best(k)*T1(k); zs_j = best(k-1)*T1(k-1);
    c_str = (zs_k - zs_j)/(T1(k) - T1(k-1));
    c_dix = sqrt((best(k)^2*T1(k) - best(k-1)^2*T1(k-1))/(T1(k) - T1(k-1)));
    fprintf('  %4.0f-%-4.0f      %10.1f          %10.1f\n', zt(k-1)*1e3, zt(k)*1e3, c_str, c_dix);
end
k40 = find(zt == 40e-3); w1 = (d_fat/c_fat)/T1(k40); w2 = 1 - w1;
fprintf('Prediction V_rms - V_h ~ w1 w2 (c_liv - c_fat)^2/(2 V_h) at 40 mm: %.2f m/s (exact %.2f)\n', ...
    w1*w2*(c_liv - c_fat)^2/(2*Vh(k40)), Vr(k40) - Vh(k40));

figure('Position', [50 50 1000 420]);
subplot(1, 2, 1);
plot(chat, CF(zt == 15e-3, :)/max(CF(zt == 15e-3, :)), 'b', chat, CF(zt == 30e-3, :)/max(CF(zt == 30e-3, :)), 'r', ...
     chat, CF(zt == 50e-3, :)/max(CF(zt == 50e-3, :)), 'k', 'LineWidth', 1.4);
xlabel('assumed sound speed c_{hat} (m/s)'); ylabel('CF / max'); legend('z = 15 mm', 'z = 30 mm', 'z = 50 mm', 'Location', 'southwest');
title('Coherence factor versus assumed speed'); grid on;
subplot(1, 2, 2);
plot(best, zt*1e3, 'ko', Vh, zt*1e3, 'b-', Vr, zt*1e3, 'r--', 'LineWidth', 1.4); set(gca, 'YDir', 'reverse');
xlabel('speed (m/s)'); ylabel('depth (mm)'); legend('CF-optimal c_{hat}', 'V_h = z/T', 'V_{rms} (Dix)', 'Location', 'southwest');
title('The focus-optimal speed is an average'); grid on;
