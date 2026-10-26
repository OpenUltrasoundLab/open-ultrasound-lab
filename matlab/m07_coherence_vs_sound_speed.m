%% m07_coherence_vs_sound_speed.m  --  toy coherence-based sound-speed estimator
% Purpose : simulate narrowband focused channel data from diffuse scatterers in a medium
%           of sound speed c_true, then re-focus (transmit and receive delays) with an
%           assumed speed c_hat swept over 1440-1640 m/s. The coherence factor and the
%           short-lag spatial coherence (SLSC) both peak where c_hat = c_true, because any
%           residual (quadratic) focusing phase can only lower spatial coherence.
%           Repeat with a near-field phase screen and with noise to see the peak drop,
%           flatten and (for the screen) shift.
% Change  : c_true, z (focal depth), N (aperture), sig_phi, SNR_dB, M (realizations).
% Expect  : clean argmax ~ c_true for both metrics; SLSC loses ~4-5 % and CF ~30-36 % at
%           |c_hat - c_true| = 40 m/s (F/1.6); the screen lowers and flattens the peak and
%           can move it by tens of m/s (here SLSC +25 m/s): aberration biases the estimate.
% Model   : narrowband phasors (single frequency f0), scatterers on a fine grid in the
%           plane at the true depth z; the beamformer believes the echo comes from
%           z_hat = z * c_hat / c_true and focuses there with speed c_hat.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(2);

%% 1. Geometry and data
c_true = 1510;  f0 = 5e6;  k = 2*pi*f0/c_true;       % true speed (m/s), Hz, wavenumber (rad/m)
N = 64;  pitch = 0.3e-3;  z = 30e-3;                  % elements, m, true focal depth (m)
xe = ((0:N-1) - (N-1)/2) * pitch;                     % element positions (m)
xi = -6e-3 : 20e-6 : 6e-3;  Nxi = numel(xi);          % scatterer grid in the focal plane (m)
M  = 300;                                             % speckle realizations (snapshots)
H  = exp(-1i*k*sqrt((xe.' - xi).^2 + z^2));           % true one-way propagation, N x Nxi
A  = (randn(Nxi, M) + 1i*randn(Nxi, M)) / sqrt(2);   % reflectivity, same draws for every c_hat

sig_phi = 0.6;  ell = 5;                              % phase screen: rms (rad), corr. length (el.)
w = conv(randn(1, N + 40), exp(-(-20:20).^2 / ell^2), 'same');  phi = w(21:20+N);
phi = phi - mean(phi);  phi = sig_phi * phi / sqrt(mean(phi.^2));
SNR_dB = 10;                                          % channel SNR for the noisy case

Mlag = 6;
lagavg = @(C, mmax) arrayfun(@(m) mean(real(diag(C, m))), 1:mmax);
cvec = 1440:5:1640;  nc = numel(cvec);
CF = zeros(3, nc);  SL = zeros(3, nc);                % rows: clean, screen, noise
Wn = (randn(N, M) + 1i*randn(N, M)) / sqrt(2);       % fixed noise draw

%% 2. Sweep the assumed sound speed
for ic = 1:nc
    ch   = cvec(ic);
    zh   = z * ch / c_true;                            % depth the beamformer assigns to the echo
    foc  = exp(1i*2*pi*f0*sqrt(xe.'.^2 + zh^2)/ch);    % focusing phase (delay) per element, N x 1
    for cs = 1:3
        if cs == 2, Hc = H .* exp(1i*phi.'); else, Hc = H; end
        Pt = (foc.' * Hc);                             % transmit field at the scatterers, 1 x Nxi
        S  = (foc .* Hc) .* Pt;                        % receive focusing x propagation x transmit
        S  = S * A;                                    % N x M focused channel samples
        if cs == 3
            Ps = mean(abs(S(:)).^2);  S = S + sqrt(Ps / 10^(SNR_dB/10)) * Wn;
        end
        R  = (S * S') / M;  d = sqrt(real(diag(R)));  C = R ./ (d * d.');
        coh = lagavg(C, Mlag);
        SL(cs, ic) = sum(coh);
        CF(cs, ic) = mean(abs(sum(S, 1)).^2) / mean(N * sum(abs(S).^2, 1));
    end
end

%% 3. Peak picking with a parabolic fit around the maximum
pfit = @(ym, y0, yp) 0.5 * (ym - yp) / (ym - 2*y0 + yp);   % sub-step offset of the vertex
names = {'clean', sprintf('screen %.1f rad', sig_phi), sprintf('SNR %d dB', SNR_dB)};
ideal = sum(1 - (1:Mlag)/N);
i0 = find(cvec == c_true);  ip = find(cvec == c_true + 40);  im = find(cvec == c_true - 40);
for cs = 1:3
    pk = zeros(1, 2);  Y = {CF(cs,:), SL(cs,:)};
    for q = 1:2
        y = Y{q};  ib = min(max(find(y == max(y), 1), 2), nc-1);
        pk(q) = cvec(ib) + 5 * pfit(y(ib-1), y(ib), y(ib+1));
    end
    fprintf('%-16s argmax CF %6.1f m/s  SLSC %6.1f m/s | peak CF %.3f SLSC/ideal %.3f | at +-40 m/s: CF %.0f%%/%.0f%%, SLSC %.1f%%/%.1f%% of peak\n', ...
        names{cs}, pk(1), pk(2), max(CF(cs,:)), max(SL(cs,:))/ideal, ...
        100*CF(cs,im)/max(CF(cs,:)), 100*CF(cs,ip)/max(CF(cs,:)), 100*SL(cs,im)/max(SL(cs,:)), 100*SL(cs,ip)/max(SL(cs,:)));
end
lam = c_true / f0;  D = N * pitch;
fprintf('edge defocus phase for 40 m/s, paraxial estimate pi D^2/(2 lambda z) * dc/c = %.2f rad\n', ...
    pi * D^2 / (2 * lam * z) * 40 / c_true);

%% 4. Figures
figure('Position', [100 100 1050 400]);
subplot(1,2,1);
plot(cvec, CF(1,:), 'b-', 'LineWidth', 1.6); hold on;
plot(cvec, CF(2,:), 'Color', [0 0.6 0.4], 'LineWidth', 1.6);
plot(cvec, CF(3,:), 'r-', 'LineWidth', 1.6);
plot([c_true c_true], [0 1], 'k--');
grid on; ylim([0 0.75]); xlim([1440 1640]); xlabel('assumed sound speed (m/s)'); ylabel('coherence factor (ratio of means)');
legend(names{1}, names{2}, names{3}, 'true c', 'Location', 'south'); title('CF vs assumed sound speed');
subplot(1,2,2);
plot(cvec, SL(1,:)/ideal, 'b-', 'LineWidth', 1.6); hold on;
plot(cvec, SL(2,:)/ideal, 'Color', [0 0.6 0.4], 'LineWidth', 1.6);
plot(cvec, SL(3,:)/ideal, 'r-', 'LineWidth', 1.6);
plot([c_true c_true], [0 1.05], 'k--');
grid on; ylim([0 1.05]); xlim([1440 1640]); xlabel('assumed sound speed (m/s)'); ylabel('SLSC (lags 1-6) / ideal');
title('SLSC vs assumed sound speed');
