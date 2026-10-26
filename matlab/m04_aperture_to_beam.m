%% m04_aperture_to_beam.m  -- Module 04: the beam pattern is the Fourier transform of the aperture
% Purpose : compute the far-field (or focal-plane) pattern of a linear array as the DTFT of
%           its apodization, compare rect vs Hann, one-way vs two-way (Tx = Rx), show that a
%           zero-padded FFT of the weights reproduces the direct far-field sum, and move
%           grating lobes with the pitch (d = lambda/2, lambda, 2 lambda).
% Change  : f0, N (elements), the pitch list, the steering angle th0_deg.
% Expect  : first sidelobe  rect -13.3 dB one-way / -26.5 dB two-way,
%                           Hann -31.5 dB one-way / -63 dB two-way;
%           -6 dB width (units of lambda/D in sin(theta)): rect 1.21 / 0.89, Hann 2.00 / 1.44;
%           FFT vs direct sum agree to ~1e-14; grating lobes at sin(th_g) = sin(th0) + m lambda/d
%           (d = lambda, 20 deg steering: -41.1 deg); element factor (w = 0.9 d) lowers that
%           grating lobe by ~4.4 dB one-way; a 60 % band fills the rect two-way first null
%           (-24 dB) and lowers the level at the narrowband sidelobe position to -28.8 dB.
% Model   : sections 1-3 monochromatic, point-like elements (no element directivity), far field
%           or focal plane; section 4 adds a rigid-baffle element factor sinc(w sin(theta)/lambda);
%           section 5 a Gaussian two-way spectrum (time-integrated energy pattern).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

%% Parameters
c   = 1.540;            % mm/us
f0  = 5;                % MHz
lam = c/f0;             % wavelength (mm) = 0.308
N   = 64;               % elements
d   = lam/2;            % pitch for the apodization study (mm)
n   = 0:N-1;
xn  = (n - (N-1)/2)*d;  % element positions (mm), centered
D   = N*d;              % aperture length (mm)
wR  = ones(1, N);
wH  = 0.5 - 0.5*cos(2*pi*(n + 0.5)/N);       % symmetric Hann sampled at element centers
AF  = @(w, x, u) abs( exp(1j*2*pi/lam * u(:) * x) * w(:) ).' / sum(w);   % |sum w e^{jk x u}|, u = sin(th) - sin(th0)

%% 1. First sidelobe and -6 dB width: rect vs Hann, one-way vs two-way
uf = linspace(0, 6*lam/D, 6001);              % fine grid from the main-lobe peak outward
fprintf('--- apodization (N = %d, d = lambda/2, D = %.2f mm) ---\n', N, D);
fprintf('%-10s %-8s %12s %16s %14s\n', 'aperture', 'pattern', '1st SL (dB)', '-6 dB width', 'first null');
labels = {'rect', 'Hann'}; Ws = {wR, wH};
theorySL = [-13.26 -26.52; -31.47 -62.93]; theoryW6 = [1.207 0.886; 2.000 1.441];
for ia = 1:2
    a1 = AF(Ws{ia}, xn, uf);
    for way = 1:2
        a = a1.^way; dBv = 20*log10(a + 1e-300);
        k = 2; while ~(dBv(k) < dBv(k-1) && dBv(k) <= dBv(k+1)), k = k + 1; end      % first null
        j = k; while ~(dBv(j) > dBv(j-1) && dBv(j) >= dBv(j+1)), j = j + 1; end      % first sidelobe
        i6 = find(dBv < -6.0206, 1);
        u6 = uf(i6-1) + (-6.0206 - dBv(i6-1))*(uf(i6) - uf(i6-1))/(dBv(i6) - dBv(i6-1));
        fprintf('%-10s %-8s %7.2f (%6.2f) %9.3f (%5.3f) lam/D %8.3f lam/D\n', labels{ia}, char('one-way' * (way == 1) + 'two-way' * (way == 2)), ...
                dBv(j), theorySL(ia, way), 2*u6*D/lam, theoryW6(ia, way), uf(k)*D/lam);
    end
end
fprintf('(lateral -6 dB width in the focal plane = width in lam/D times lambda F#: rect two-way 0.886 lambda F#)\n');
Fn = 2;                                       % example F-number
fprintf('at F# = %d, %.0f MHz (lambda F# = %.3f mm): rect -6 dB width %.3f mm one-way, %.3f mm two-way\n', ...
        Fn, f0, lam*Fn, 1.207*lam*Fn, 0.886*lam*Fn);
fprintf('Hann/rect -6 dB width ratio: %.2f one-way, %.2f two-way\n', 2.000/1.207, 1.441/0.886);

%% 2. Beam = FFT of the aperture: zero-padded FFT vs direct far-field sum
L = 4096;
Wf = abs(fft(wH, L)) / sum(wH);                % DTFT samples at psi = -2 pi k/L
psi = -2*pi*(0:L-1)/L; psi = mod(psi + pi, 2*pi) - pi;     % wrap to [-pi, pi)
uk = psi*lam/(2*pi*d);                         % psi = k d u  ->  u = sin(theta)
vis = abs(uk) <= 1;
direct = AF(wH, xn, uk(vis));
fprintf('\n--- FFT of the weights vs direct far-field sum (Hann, L = %d) ---\n', L);
fprintf('max | |FFT| - |direct sum| | over %d visible-region samples = %.2e\n', nnz(vis), max(abs(Wf(vis) - direct)));
fprintf('sin(theta) per FFT bin = lambda/(L d) = %.5f\n', lam/(L*d));

%% 3. Grating lobes vs pitch
th0_deg = 20;
pitches = [0.5 1 2]*lam;
uu = linspace(-1, 1, 20001); thd = asind(uu);
fprintf('\n--- grating lobes, N = %d, steering 0 and %d deg ---\n', N, th0_deg);
figure('Name', 'grating lobes', 'Position', [80 80 760 640]);
for ip = 1:3
    dp = pitches(ip); xp = (n - (N-1)/2)*dp;
    for st = [0 th0_deg]
        s0 = sind(st);
        a = AF(wR, xp, uu - s0); adB = 20*log10(a + 1e-300);
        pk = find([adB(1) > adB(2), adB(2:end-1) > adB(1:end-2) & adB(2:end-1) >= adB(3:end), adB(end) > adB(end-1)] & adB > -1.5);
        meas = thd(pk); meas = meas(abs(meas - st) > 0.5);                          % drop the main lobe
        m = -4:4; m = m(m ~= 0); sg = s0 + m*lam/dp; pred = asind(sg(abs(sg) <= 1 + 1e-12));
        sp = sprintf(' %.1f', sort(pred)); if isempty(pred), sp = ' none'; end
        sm = sprintf(' %.1f', sort(meas)); if isempty(meas), sm = ' none'; end
        fprintf('d = %.2f lambda (%.3f mm), steer %2d deg: predicted [%s ] deg, measured [%s ] deg\n', dp/lam, dp, st, sp, sm);
        if st == th0_deg
            subplot(3, 1, ip); plot(thd, max(-60, adB), 'b'); hold on;
            for g = pred, plot([g g], [-60 3], 'r--'); end
            plot([st st], [-60 3], 'k:'); hold off; ylim([-60 3]); xlim([-90 90]); grid on;
            xlabel('angle \theta (deg)'); ylabel('|AF| (dB)');
            title(sprintf('d = %.1f \\lambda, steered %d deg: red = predicted grating lobes', dp/lam, st));
        end
    end
end
fprintf('d <= lambda/2: no grating lobe for any steering; d <= lambda/(1+sin(th_max)) for |th0| <= th_max\n');
fprintf('for th_max = %d deg: d <= %.3f mm (%.3f lambda)\n', th0_deg, lam/(1 + sind(th0_deg)), 1/(1 + sind(th0_deg)));

%% 4. Element directivity attenuates grating lobes (d = lambda, steered 20 deg)
dp = lam; xp = (n - (N-1)/2)*dp; s0 = sind(th0_deg);
we = 0.9*dp;                                       % element width (mm): 10 % kerf
sincf = @(v) (sin(pi*v) + (v == 0)) ./ (pi*v + (v == 0));
EF = @(uabs) abs(sincf(we*uabs/lam));
thg = asind(s0 - lam/dp);                          % predicted grating lobe (deg)
aM = AF(wR, xp, 0)        * EF(s0);                % main lobe (u = 0 at theta = th0)
aG = AF(wR, xp, sind(thg) - s0) * EF(sind(thg));   % grating lobe
fprintf('\n--- element directivity (d = lambda, w = 0.9 d, steer %d deg) ---\n', th0_deg);
fprintf('grating lobe at %.1f deg; element factor there %.2f dB, at the main lobe %.2f dB\n', thg, 20*log10(EF(sind(thg))), 20*log10(EF(s0)));
fprintf('grating lobe re main lobe: %.2f dB one-way, %.2f dB two-way (0 dB without the element factor)\n', 20*log10(aG/aM), 40*log10(aG/aM));

%% 5. Broadband: a 60 % two-way band smears sidelobes and fills nulls (energy pattern)
FBW = 0.60; sf = FBW/(2*sqrt(2*log(2)));           % Gaussian amplitude-spectrum std (units of f0)
fr  = linspace(1 - 4*sf, 1 + 4*sf, 41);            % frequencies (units of f0)
Pw  = exp(-(fr - 1).^2/sf^2);                      % |P(f)|^2 weights
ue  = [0 1 1.43 2 2.46]*lam/D;                     % sin(theta): peak, rect null, rect SL, Hann null, Hann SL
AFl = @(w, x, u, lm) abs( exp(1j*2*pi/lm * u(:) * x) * w(:) ).' / sum(w);
fprintf('\n--- broadband two-way energy pattern, FBW = %.0f %% (levels in dB re peak) ---\n', 100*FBW);
fprintf('%-6s %-11s %8s %8s %8s %8s\n', 'apod', 'model', 'u=1', 'u=1.43', 'u=2', 'u=2.46');
ub2 = linspace(-8*lam/D, 8*lam/D, 4001); Ebb = zeros(2, numel(ub2));
for ia = 1:2
    Eu = zeros(1, numel(ue)); Eb = zeros(1, numel(ub2));
    for kf = 1:numel(fr)
        Eu = Eu + Pw(kf)*AFl(Ws{ia}, xn, ue,  lam/fr(kf)).^4;
        Eb = Eb + Pw(kf)*AFl(Ws{ia}, xn, ub2, lam/fr(kf)).^4;
    end
    Ebb(ia, :) = 10*log10(Eb/max(Eb));
    nb = 40*log10(AF(Ws{ia}, xn, ue) + 1e-300);
    fprintf('%-6s %-11s %8.1f %8.1f %8.1f %8.1f\n', labels{ia}, 'narrowband', max(nb(2:end), -99));
    fprintf('%-6s %-11s %8.1f %8.1f %8.1f %8.1f\n', labels{ia}, 'broadband', 10*log10(Eu(2:end)/Eu(1)));
end
fprintf('(u in units of lambda0/D; -99 marks an exact null)\n');

%% Figures: apodization and FFT-vs-direct
figure('Name', 'aperture to beam', 'Position', [80 80 760 640]);
subplot(2, 2, 1); plot(xn, wR, 'bo', xn, wH, 'r.-', 'MarkerSize', 8); ylim([0 1.1]);
xlabel('element position (mm)'); ylabel('weight'); title('aperture: rect (blue), Hann (red)');
ub = linspace(-8*lam/D, 8*lam/D, 4001);
subplot(2, 2, 2);
plot(ub*D/lam, 20*log10(AF(wR, xn, ub)), 'b', ub*D/lam, 20*log10(AF(wH, xn, ub)), 'r'); ylim([-80 3]); grid on;
xlabel('sin\theta in units of \lambda/D'); ylabel('dB'); title('one-way: -13.3 vs -31.5 dB');
subplot(2, 2, 3);
plot(ub*D/lam, 40*log10(AF(wR, xn, ub)), 'b', ub*D/lam, 40*log10(AF(wH, xn, ub)), 'r'); hold on;
plot(ub2*D/lam, Ebb(1, :), 'b--', ub2*D/lam, Ebb(2, :), 'r--', 'LineWidth', 1.2); hold off; ylim([-100 3]); grid on;
xlabel('sin\theta in units of \lambda/D'); ylabel('dB'); title('two-way: narrowband (solid), 60 % band (dashed)');
subplot(2, 2, 4);
ud = linspace(-1, 1, 4001);
plot(ud, 20*log10(AF(wH, xn, ud)), 'k'); hold on;
iv = find(vis); iv = iv(1:4:end);
plot(uk(iv), 20*log10(Wf(iv) + 1e-300), 'ro', 'MarkerSize', 3); hold off; ylim([-100 3]); xlim([-0.3 0.3]); grid on;
xlabel('sin\theta'); ylabel('dB'); title('Hann: direct sum (black) vs FFT bins (red)');
