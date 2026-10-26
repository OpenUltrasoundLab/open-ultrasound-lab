%% m14_autofocus_sound_speed.m  --  sound speed by gradient ascent on a focus metric
% Purpose : simulate plane-wave IQ channel data (64 elements, 5 MHz) from a point
%           target in speckle in a medium with c_true = 1480 m/s, beamform at an
%           assumed speed ch, and estimate c by gradient ascent on four metrics:
%             energy : sum of |b|^2 over an ROI fixed in two-way time (coherently compounded)
%             sharp  : P*sum(I.^2)/sum(I)^2 over the same ROI
%             envz   : |b|^2 at a pixel fixed in depth (x = 0, z = 30 mm)
%             rfz    : real part of b (the RF value) at that fixed pixel
%           The derivative dM/dch is computed exactly by the chain rule (forward mode):
%           d/dch [IQ(tau) e^{j w0 tau}] = (IQ'(tau) + j w0 IQ(tau)) e^{j w0 tau} dtau/dch.
% Change  : c_true, z, angles_deg (try a single 0 deg plane wave), starts, the scene.
% Expect  : 5 angles: energy peaks at 1484-1486 m/s (bias about +5 m/s) and every start
%           in 1380-1720 converges there; rfz oscillates with period ch^2/(2 z f0)
%           (7.3 m/s here) and every start stops on its nearest local maximum (cycle
%           skipping). 1 angle: energy is nearly flat and peaks far from c_true.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8. Runtime about 10 s.
rng(1);

%% 1. Geometry, pulse and scatterers
c_true = 1480;            % true sound speed (m/s)
z      = 30e-3;           % target depth (m)
f0 = 5e6; w0 = 2*pi*f0;   % centre frequency (Hz)
sig = 0.125e-6;           % Gaussian envelope std (s): -6 dB bandwidth 3 MHz (60 %)
fs = 40e6;                % IQ sampling rate (Hz)
N = 64; pitch = 0.3e-3;   % elements, pitch (m): aperture 19.2 mm, F/1.56 at 30 mm
xe = ((0:N-1) - (N-1)/2) * pitch;                   % 1 x N element positions (m)
K = 500;                                            % speckle scatterers
sx = [0; (rand(K,1) - 0.5) * 8e-3];                 % lateral (m), first one is the point
sz = [z; z + (rand(K,1) - 0.5) * 3e-3];             % depth (m)
amp = [3; 0.25 * (randn(K,1) + 1j*randn(K,1)) / sqrt(2)];
angle_sets = {0, (-8:4:8)};                         % plane-wave angles (deg)

for set_id = 1:2
th = angle_sets{set_id} * pi/180;  A = numel(th);

%% 2. IQ channel data: IQ(t) = sum_k a_k g(t - t_k) e^{-j w0 t_k}
tk = zeros(K+1, N, A);
for a = 1:A
    tk(:,:,a) = (sz*cos(th(a)) + sx*sin(th(a)))/c_true + sqrt((sx - xe).^2 + sz.^2)/c_true;
end
t0 = min(tk(:)) - 1.5e-6;  T = ceil((max(tk(:)) + 1.5e-6 - t0) * fs);
tg = t0 + (0:T-1).' / fs;                            % T x 1 fast time (s)
IQ = zeros(T, N, A);
for a = 1:A
    for e = 1:N
        tt = tk(:, e, a).';                          % 1 x (K+1)
        IQ(:, e, a) = exp(-(tg - tt).^2 / (2*sig^2)) * (amp .* exp(-1j*w0*tt.'));
    end
end

%% 3. Beamformer with exact derivative (anonymous functions only)
% ROI in time coordinates: pixel (x, t), depth at assumed speed = ch*t/2
nx = 21; nt = 25;
[PX, PT] = meshgrid(((0:nx-1) - (nx-1)/2) * 0.15e-3, 2*z/c_true + ((0:nt-1) - (nt-1)/2) * 50e-9);
Px = PX(:); Pt = PT(:); P = numel(Px);
U   = @(tau) (tau - t0) * fs;
I0  = @(tau) min(max(floor(U(tau)), 0), T-2);
% core: q0 = IQ at the sample before tau, dq = next minus q0, fr = fraction, E = e^{j w0 tau}
core = @(q0, dq, fr, E, dt) [sum((q0 + dq.*fr) .* E, 2), sum((dq*fs + 1j*w0*(q0 + dq.*fr)) .* E .* dt, 2)];
bfi  = @(k, fr, tau, dt) core(IQ(k), IQ(k+1) - IQ(k), fr, exp(1j*w0*tau), dt);   % k = linear index
bfa  = @(tau, dt, a) bfi(I0(tau) + 1 + (0:N-1)*T + (a-1)*T*N, U(tau) - I0(tau), tau, dt);
% time-ROI delays and d tau / d ch
tauT  = @(ch, a) Pt/2*cos(th(a)) + Px*sin(th(a))/ch + sqrt((Px - xe).^2/ch^2 + (Pt/2).^2);
dtauT = @(ch, a) -Px*sin(th(a))/ch^2 - ((Px - xe).^2/ch^3) ./ sqrt((Px - xe).^2/ch^2 + (Pt/2).^2);
% fixed-depth pixel (x = 0, z)
tauZ  = @(ch, a) (z*cos(th(a)) + sqrt(xe.^2 + z^2)) / ch;
BT = @(ch) cell2mat(arrayfun(@(a) bfa(tauT(ch,a), dtauT(ch,a), a), 1:A, 'UniformOutput', false));
BZ = @(ch) cell2mat(arrayfun(@(a) bfa(tauZ(ch,a), -tauZ(ch,a)/ch, a), 1:A, 'UniformOutput', false));
bsum = @(Bm) [sum(Bm(:,1:2:end), 2), sum(Bm(:,2:2:end), 2)];       % compound angles: [b, db]
% metrics: each returns [M, dM/dch]
mE = @(bd) [sum(abs(bd(:,1)).^2), sum(2*real(conj(bd(:,1)).*bd(:,2)))];
mS = @(bd) P * [sum(abs(bd(:,1)).^4) / sum(abs(bd(:,1)).^2)^2, ...
     (sum(4*abs(bd(:,1)).^2 .* real(conj(bd(:,1)).*bd(:,2))) * sum(abs(bd(:,1)).^2) ...
      - 2*sum(abs(bd(:,1)).^4) * sum(2*real(conj(bd(:,1)).*bd(:,2)))) / sum(abs(bd(:,1)).^2)^3];
metric = {@(ch) mE(bsum(BT(ch))), @(ch) mS(bsum(BT(ch))), ...
          @(ch) mE(bsum(BZ(ch))), @(ch) real(bsum(BZ(ch)))};
names = {'energy', 'sharp', 'envz', 'rfz'};
S0 = [40 40 8 1.5];       % initial trial step of the line search (m/s)
cmin = 1350; cmax = 1750; % box constraint: every step is projected onto [cmin, cmax]

%% 4. Derivative check at 1540 m/s (exact chain rule vs central difference)
if set_id == 2
    for m = 1:4
        r = metric{m}(1540);  h = 1e-3;
        fd = (metric{m}(1540+h) - metric{m}(1540-h)) / (2*h);
        fprintf('dM/dc check %-6s  chain rule / central difference = %.6f\n', names{m}, r(2)/fd(1));
    end
end

%% 5. Sweep and gradient ascent from several starts
fprintf('\n%d plane wave(s): c_true = %d m/s, z = %d mm\n', A, c_true, round(z*1e3));
fprintf('metric | argmax (m/s) | gradient ascent: start -> end (iterations)\n');
starts = [1380 1440 1540 1620 1720];
cs_roi = 1350:5:1750;  cs_pix = 1350:1:1750;
curves = cell(1,4);  ends = zeros(4, numel(starts));
for m = 1:4
    if m <= 2, cs = cs_roi; else, cs = cs_pix; end
    M = zeros(size(cs));  for i = 1:numel(cs), r = metric{m}(cs(i)); M(i) = r(1); end
    curves{m} = [cs; M];  [~, im] = max(M);
    s_out = '';
    for j = 1:numel(starts)
        c = starts(j);  r = metric{m}(c);  n = 0;
        while n < 60
            g = r(2);  if g == 0, break; end
            s = S0(m);  ok = false;
            while s > 0.05
                cn = min(max(c + sign(g)*s, cmin), cmax);  rn = metric{m}(cn);   % projected step
                if cn ~= c && rn(1) >= r(1) + 1e-4*abs(cn-c)*abs(g), c = cn; r = rn; ok = true; break; end
                s = s/2;
            end
            if ~ok, break; end
            n = n + 1;
        end
        ends(m, j) = c;
        s_out = [s_out sprintf('%d->%.1f(%d)  ', starts(j), c, n)];
    end
    fprintf('%-6s | %7.1f      | %s\n', names{m}, cs(im), s_out);
end
fprintf('RF-at-fixed-depth period predicted ch^2/(2 z f0) = %.1f m/s at ch = c_true\n', c_true^2/(2*z*f0));
if set_id == 1, curves1 = curves; ends1 = ends; else, curves5 = curves; ends5 = ends; end
end

%% 6. Figures
figure('Position', [80 80 1100 460]);
cols = {[0 0.3 0.8], [0.8 0.2 0.5], [0.85 0.55 0], [0 0.55 0.4]};
for set_id = 1:2
    subplot(1,2,set_id);
    if set_id == 1, cv = curves1; en = ends1; ttl = '1 plane wave (0 deg)'; else, cv = curves5; en = ends5; ttl = '5 plane waves (-8:4:8 deg)'; end
    h = zeros(1,4);
    for m = 1:4
        y = cv{m}(2,:);  y = y / max(abs(y));
        h(m) = plot(cv{m}(1,:), y, 'Color', cols{m}, 'LineWidth', 1.3); hold on;
        yi = interp1(cv{m}(1,:), y, en(m,:));
        plot(en(m,:), yi, 'o', 'Color', cols{m}, 'MarkerSize', 7, 'LineWidth', 1.5);
    end
    plot([c_true c_true], [-1 1.05], 'k--');
    xlim([1350 1750]); ylim([-1 1.05]); grid on;
    xlabel('assumed sound speed c_h (m/s)'); ylabel('metric / max |metric|');
    title([ttl ': circles = gradient-ascent end points']);
    legend(h, {'energy (time ROI)', 'sharpness (time ROI)', '|b|^2 at fixed depth', 'Re b at fixed depth'}, 'Location', 'southwest');
end
