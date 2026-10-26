%% m14_cycle_skipping.m  --  why RF misfits trap gradient descent, and how to widen the basin
% Purpose : delay estimation by minimizing L(tau) = sum (d(t) - s(t - tau))^2 for a 5 MHz
%           Gaussian pulse (60 % bandwidth). Compares three misfits computed numerically
%           on sampled signals: RF, envelope (|analytic signal|, FFT Hilbert), and RF after
%           a Gaussian low-pass at 1.5 or 1.0 MHz. Lists the local minima, runs gradient
%           descent from 121 starting delays, and shows that a multiscale schedule (envelope
%           first, then RF) converges from starts the RF misfit alone cannot handle. Ends with
%           the travel-time form of the criterion used in full-waveform inversion.
% Change  : f0, sg, the low-pass corner fl, the starting delays.
% Expect  : RF local minima at 0, +-194, +-388 ns (just inside multiples of T0 = 200 ns);
%           RF basin |tau0| < 97 ns; envelope converges from every start (no local minima);
%           low-pass at 1.5 MHz: basin about +-160 ns at -29 dB of the pulse energy;
%           envelope first, then RF: converges from every start with RF curvature (32x
%           the envelope curvature at the minimum).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8. Runtime about 5 s.
rng(1);

%% 1. Signals
f0 = 5e6; sg = 0.125e-6;           % centre frequency (Hz), envelope std (s)
fs = 1e9; t = (-3e-6:1/fs:3e-6).';  % fine grid so delays are continuous (1 ns)
nT = numel(t);  nf = 2^nextpow2(2*nT);
fax = [0:nf/2-1, -nf/2:-1].' * fs / nf;              % FFT frequency axis (Hz)
pulse = @(tau) exp(-(t-tau).^2/(2*sg^2)) .* cos(2*pi*f0*(t-tau));
fl = 1.5e6;                                           % Gaussian low-pass corner (Hz)
lowp  = @(x) real(ifft(fft(x, nf) .* exp(-fax.^2/(2*fl^2))));
H     = @(x) ifft(fft(x, nf) .* (1 + sign(fax)));     % analytic signal (x + j Hilbert x)
crop  = @(x) x(1:nT);
envl  = @(x) abs(crop(H(x)));
d = pulse(0);                                          % data: true delay 0
fprintf('low-pass at %.1f MHz keeps %.1f dB of the pulse energy\n', fl/1e6, 10*log10(sum(crop(lowp(d)).^2)/sum(d.^2)));

%% 2. Misfit curves, computed numerically from the sampled signals
% L(tau) = sum (a - b_tau)^2 = sum a^2 + sum b^2 - 2 R_ab(tau): all shifts at once by FFT
% correlation (a shift commutes with the low-pass filter and with taking the envelope).
xc = @(a, b) real(ifft(fft(a, nf) .* conj(fft(b, nf))));
lagk = -600:600;  taus = lagk * 1e-9;                  % 1 ns lags
mis = @(a, b) sum(a.^2) + sum(b.^2) - 2 * subsref(xc(b, a), struct('type', '()', 'subs', {{mod(lagk, nf) + 1}})).';
dE = envl(d);  dL = crop(lowp(d));
fl2 = 1.0e6;  lowp2 = @(x) real(ifft(fft(x, nf) .* exp(-fax.^2/(2*fl2^2))));  dL2 = crop(lowp2(d));
L = [mis(d, d); mis(dE, dE); mis(dL, dL); mis(dL2, dL2)];
L = L ./ (2*[sum(d.^2); sum(dE.^2); sum(dL.^2); sum(dL2.^2)]);   % normalize by 2 x energy: L = 1 - R(tau)/R(0)
fprintf('low-pass at %.1f MHz keeps %.1f dB of the pulse energy\n', fl2/1e6, 10*log10(sum(dL2.^2)/sum(d.^2)));
names = {'RF', 'envelope', 'LP 1.5 MHz', 'LP 1.0 MHz'};
for m = 1:4
    loc = find(L(m,2:end-1) < L(m,1:end-2) & L(m,2:end-1) < L(m,3:end)) + 1;
    fprintf('%-11s local minima (ns): %s\n', names{m}, sprintf('%g ', round(taus(loc)*1e9)));
end
mx = find(L(1,2:end-1) > L(1,1:end-2) & L(1,2:end-1) > L(1,3:end)) + 1;
fprintf('RF local maxima nearest 0 (basin edges): %s ns; T0/2 = %.0f ns\n', sprintf('%g ', round(taus(mx(abs(taus(mx)) < 150e-9))*1e9)), 0.5e9/f0);

%% 3. Gradient descent from many starts on the tabulated misfits
% dL/dtau by central differences on the 1 ns table; step = 1/L''(0) (Newton step at the minimum)
dLdt = zeros(size(L));  dLdt(:,2:end-1) = (L(:,3:end) - L(:,1:end-2)) / 2e-9;
curv = (L(:,602) - 2*L(:,601) + L(:,600)) / 1e-18;    % L'' at tau = 0
fprintf('curvature at the minimum, RF / envelope = %.1f (theory 1 + 2 (2 pi f0 sg)^2 = %.1f)\n', curv(1)/curv(2), 1 + 2*(2*pi*f0*sg)^2);
starts = (-600:10:600) * 1e-9;
plan = {1, 2, 3, 4, [2 1]};                           % misfits used in sequence
lbl = {'RF', 'envelope', 'LP 1.5 MHz', 'LP 1.0 MHz', 'envelope -> RF'};
final = zeros(numel(plan), numel(starts));
for m = 1:numel(plan)
    tau = starts;
    for q = plan{m}
        for it = 1:400
            g = interp1(taus, dLdt(q,:), min(max(tau, taus(2)), taus(end-1)));
            tau = min(max(tau - g / curv(q), taus(1)), taus(end));   % stay on the table
        end
    end
    final(m, :) = tau;
end
ok = abs(final) < 2e-9;
for m = 1:numel(plan)
    k0 = find(starts == 0);  lo = k0;  hi = k0;
    while lo > 1 && ok(m, lo-1), lo = lo - 1; end
    while hi < numel(starts) && ok(m, hi+1), hi = hi + 1; end
    fprintf('%-15s converges to 0 from %3d of %d starts; contiguous basin [%4.0f, %4.0f] ns\n', ...
        lbl{m}, sum(ok(m,:)), numel(starts), starts(lo)*1e9, starts(hi)*1e9);
end

%% 4. Travel-time form of the criterion (transmission FWI / tomography)
% the starting model must predict the travel time within T/2: |dt| = L |dc| / c^2 < 1/(2 f)
fprintf('\nmax tolerable starting error |dc| (m/s) for c = 1500 m/s\n   f (MHz) | L = 50 mm | L = 100 mm | L = 200 mm\n');
for f = [0.5 1 2 5]
    fprintf('   %5.1f   | %8.1f  | %9.1f  | %9.1f\n', f, 1500^2./(2*f*1e6*[0.05 0.1 0.2]));
end

%% 5. Figures
figure('Position', [100 100 1050 420]);
subplot(1,2,1);
h = plot(taus*1e9, L(1,:), 'b', taus*1e9, L(2,:), 'r', taus*1e9, L(3,:), 'g');
set(h(3), 'Color', [0 0.55 0.4]); set(h, 'LineWidth', 1.4); grid on;
xlabel('trial delay \tau (ns), true delay 0'); ylabel('misfit / (2 x signal energy)');
xlim([-600 600]); legend(h, names(1:3), 'Location', 'southeast'); title('RF misfit: a minimum every period');
subplot(1,2,2);
h2 = plot(starts*1e9, final(1,:)*1e9, 'b.-', starts*1e9, final(3,:)*1e9, 'g.-', starts*1e9, final(5,:)*1e9 + 8, 'k.-');
set(h2(2), 'Color', [0 0.55 0.4]); grid on; xlim([-600 600]); xlabel('starting delay (ns)'); ylabel('converged delay (ns)');
legend(h2, {'RF', 'LP 1.5 MHz', 'envelope -> RF (offset +8 ns)'}, 'Location', 'northwest');
title('where gradient descent ends');
