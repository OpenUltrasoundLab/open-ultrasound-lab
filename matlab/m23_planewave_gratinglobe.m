%% m23_planewave_gratinglobe.m
% Plane-wave imaging: transmit vs receive beams, lateral resolution vs angle
% span, and the grating lobe of a steered plane wave.
%
% Method (base MATLAB, no toolboxes): 128 point elements, Gaussian-envelope
% pulse (60 % band). For each plane wave the transmitted field at the target,
% h(t) = sum_e E_e * p(t - t_e - r_e/c), is synthesized on a fine time grid
% (it contains the grating-lobe wave when the pitch is too coarse). DAS of
% the analytic signal then gives, for a pixel P,
%   two-way(P) = sum_angles sum_rx w_r h(tau_tx(P) + [|r-P| - |r-T|]/c)
%   Tx-only(P) = sum_angles h(tau_tx(P))       (ideal receiver at T)
%   Rx-only(P) = sum_rx w_r p([|r-P| - |r-T|]/c) (ideal impulse at T)
% Dynamic receive aperture F# 1.5, rect, element factor sinc(w sin(theta)/lambda)
% of a 0.9d strip on transmit and receive.
%
% What to change: f0, pitch d, angle sets, target depth zt, FBW, Fnum.
% Expected output (Octave 8):
%  Part 1 (5 MHz, 0.3 mm, +-15 deg, z 30 mm): two-way -6 dB width 0.611 mm
%   for one angle (Rx-only 0.587 mm), 0.420/0.479/0.490 mm for 3/11/31 angles;
%   narrowband prediction kappa*lambda/(sin(theta_max) + 1/(2F#)) = 0.456 mm.
%  Part 3 streaks (1 angle / 11 angles): -30.2/-42.9 (5 MHz, 0.15 mm, z 15),
%   -24.3/-36.0 (5 MHz, 0.3 mm, z 15), -27.2/-35.2 (7.5 MHz, 0 deg, z 30),
%   -22.6/-35.2 dB (7.5 MHz, 15 deg, z 30).
%  Part 2: width vs span (1 deg steps) against the same formula.
%  Part 3: steered-PW grating-lobe angle, the pitch threshold
%   lambda/(1+|sin theta|), and the level of the grating-lobe streak.
% Runtime about 10 s in Octave.

clear; close all; rng(1);
c = 1.54;                 % sound speed, mm/us
Nel = 128;                % elements
FBW = 0.60;               % -6 dB fractional bandwidth of the pulse-echo pulse
Fnum = 1.5;               % receive F-number
sincf = @(x) (sin(pi*x) + (x == 0)) ./ (pi*x + (x == 0));
% kappa(rho): -6 dB width factor of sinc(2 a x) .* sinc(2 b x), computed numerically
xk = linspace(0, 6, 60001);
% lin: linear interpolation of h on its uniform grid, zero outside (q = fractional index from 0)
lin = @(h, q, L) (h(min(max(floor(q), 0), L - 2) + 1) .* (1 - (q - min(max(floor(q), 0), L - 2))) + ...
                  h(min(max(floor(q), 0), L - 2) + 2) .* (q - min(max(floor(q), 0), L - 2))) .* (q >= 0 & q <= L - 1);
w6 = @(y) 2*interp1(y(find(y < 0.5, 1) - [1 0]), xk(find(y < 0.5, 1) - [1 0]), 0.5);

%% Part 1: lateral PSF at the target depth, N angles over +-15 deg
f0 = 5; d = 0.30; zt = 30; lam = c/f0;
xs = -10:0.02:10; zs = zt*ones(size(xs));
sets = {0, linspace(-15, 15, 3), linspace(-15, 15, 11), linspace(-15, 15, 31)};
fprintf('Part 1: 5 MHz, d = %.2f mm (%.2f lambda), z = %d mm, F# %.1f\n', d, d/lam, zt, Fnum);
fprintf('  N   Tx-only -6dB   Rx-only -6dB   two-way -6dB (mm)\n');
prof = cell(1, 4);
for is = 1:numel(sets)
  ang = sets{is};
  % ---------- synthesis (same block is reused below) ----------
  xe = ((0:Nel-1) - (Nel-1)/2) * d; we = 0.9*d;
  sg = 2*sqrt(2*log(2)) / (2*pi*FBW*f0);  dt = 1/(20*f0);
  rT = sqrt(xe.^2 + zt^2); ET = sincf(we*(xe./rT)/lam);
  Tx = zeros(size(xs)); TW = zeros(size(xs));
  for a = ang
    s = sind(a); co = cosd(a);
    te = xe*s/c + rT/c; amp = ET .* sqrt(zt./rT);
    tg = (min(te) - 5*sg):dt:(max(te) + 5*sg); h = zeros(size(tg));
    for e = 1:Nel
      tt = tg - te(e); m = abs(tt) < 5*sg;
      h(m) = h(m) + amp(e) * exp(-tt(m).^2/(2*sg^2)) .* exp(2i*pi*f0*tt(m));
    end
    L = numel(tg); hi = @(t) lin(h, (t - tg(1))/dt, L);
    tau = (xs*s + zs*co)/c;
    Tx = Tx + hi(tau);
    for r = 1:Nel
      msk = abs(xe(r) - xs) <= zs/(2*Fnum);
      if ~any(msk), continue; end
      B = (sqrt((xs(msk) - xe(r)).^2 + zs(msk).^2) - rT(r))/c;
      TW(msk) = TW(msk) + ET(r) * hi(tau(msk) + B);
    end
  end
  RX = zeros(size(xs));
  for r = 1:Nel
    msk = abs(xe(r) - xs) <= zs/(2*Fnum);
    B = (sqrt((xs(msk) - xe(r)).^2 + zs(msk).^2) - rT(r))/c;
    RX(msk) = RX(msk) + ET(r) * exp(-B.^2/(2*sg^2)) .* exp(2i*pi*f0*B);
  end
  % ------------------------------------------------------------
  A = {abs(Tx)/max(abs(Tx)), abs(RX)/max(abs(RX)), abs(TW)/max(abs(TW))}; wd = zeros(1, 3);
  for k = 1:3
    y = A{k}; i0 = find(xs >= 0, 1); ir = i0 + find(y(i0:end) < 0.5, 1) - 1; il = i0 - find(y(i0:-1:1) < 0.5, 1) + 1;
    if isempty(ir) || isempty(il), wd(k) = Inf; continue; end
    xr = interp1(y([ir-1 ir]), xs([ir-1 ir]), 0.5); xl = interp1(y([il+1 il]), xs([il+1 il]), 0.5); wd(k) = xr - xl;
  end
  prof{is} = A;
  fprintf('  %2d   %8.3f       %8.3f       %8.3f\n', numel(ang), wd(1), wd(2), wd(3));
end
srx = 1/(2*Fnum); stx = sind(15);
y = abs(sincf(2*srx*xk) .* sincf(2*stx*xk)); wn = w6(y)*lam;
fprintf('  narrowband prediction, continuous +-15 deg: %.3f mm (kappa = %.3f); receive-only 1.207 lambda F# = %.3f mm\n', ...
  wn, wn/lam*(stx + srx), 1.207*lam*Fnum);

%% Part 2: width vs angular span (1 deg steps), against kappa*lambda/(sin(theta_max) + 1/(2F#))
fprintf('Part 2: two-way -6 dB width vs half-span (1 deg steps)\n  half-span  N   sim (mm)  formula (mm)  kappa\n');
for hs = [5 10 15 20]
  ang = -hs:1:hs;
  xe = ((0:Nel-1) - (Nel-1)/2) * d; we = 0.9*d; sg = 2*sqrt(2*log(2)) / (2*pi*FBW*f0); dt = 1/(20*f0);
  rT = sqrt(xe.^2 + zt^2); ET = sincf(we*(xe./rT)/lam); xs2 = -2:0.01:2; zs2 = zt*ones(size(xs2)); TW = zeros(size(xs2));
  for a = ang
    s = sind(a); co = cosd(a); te = xe*s/c + rT/c; amp = ET .* sqrt(zt./rT);
    tg = (min(te) - 5*sg):dt:(max(te) + 5*sg); h = zeros(size(tg));
    for e = 1:Nel, tt = tg - te(e); m = abs(tt) < 5*sg; h(m) = h(m) + amp(e)*exp(-tt(m).^2/(2*sg^2)).*exp(2i*pi*f0*tt(m)); end
    L = numel(tg); hi = @(t) lin(h, (t - tg(1))/dt, L);
    tau = (xs2*s + zs2*co)/c;
    for r = 1:Nel
      msk = abs(xe(r) - xs2) <= zs2/(2*Fnum); if ~any(msk), continue; end
      B = (sqrt((xs2(msk) - xe(r)).^2 + zs2(msk).^2) - rT(r))/c; TW(msk) = TW(msk) + ET(r)*hi(tau(msk) + B);
    end
  end
  y2 = abs(TW)/max(abs(TW)); i0 = find(xs2 >= 0, 1); ir = i0 + find(y2(i0:end) < 0.5, 1) - 1; il = i0 - find(y2(i0:-1:1) < 0.5, 1) + 1;
  wsim = interp1(y2([ir-1 ir]), xs2([ir-1 ir]), 0.5) - interp1(y2([il+1 il]), xs2([il+1 il]), 0.5);
  stx = sind(hs); y = abs(sincf(2*srx*xk) .* sincf(2*stx*xk)); wf = w6(y)*lam;
  fprintf('  +-%2d deg  %2d   %6.3f    %6.3f      %5.3f\n', hs, numel(ang), wsim, wf, wf/lam*(stx + srx));
end

%% Part 3: grating lobe of a steered plane wave, and its streak in the image
cases = [5 0.15 15 15; 5 0.30 15 15; 7.5 0.30 0 30; 7.5 0.30 15 30];   % f0 (MHz), d (mm), theta_tx (deg), zt (mm)
fprintf('Part 3: steered plane-wave grating lobe (single angle, then 11 angles over +-15 deg)\n');
fprintf('  f0   d/lam  th_tx  d_max/lam  th_g    Delta (mm)  streak 1 angle   11 angles (dB)\n');
img = cell(1, size(cases, 1));
for ic = 1:size(cases, 1)
  f0 = cases(ic, 1); d = cases(ic, 2); th = cases(ic, 3); zt = cases(ic, 4); lam = c/f0;
  sgl = sind(th) - lam/d; if th == 0, sgl = -lam/d; end
  if abs(sgl) <= 1, thg = asind(sgl); Dl = zt*(1 - cosd(th - thg))/cosd(thg); else thg = NaN; Dl = NaN; end
  [X, Z] = meshgrid(-10:0.15:10, (zt - 3):0.1:(zt + 9));
  lev = zeros(1, 2);
  for pass = 1:2
    if pass == 1, ang = th; else ang = linspace(-15, 15, 11); end
    xe = ((0:Nel-1) - (Nel-1)/2) * d; we = 0.9*d; sg = 2*sqrt(2*log(2)) / (2*pi*FBW*f0); dt = 1/(20*f0);
    rT = sqrt(xe.^2 + zt^2); ET = sincf(we*(xe./rT)/lam); TW = zeros(size(X));
    for a = ang
      s = sind(a); co = cosd(a); te = xe*s/c + rT/c; amp = ET .* sqrt(zt./rT);
      tg = (min(te) - 5*sg):dt:(max(te) + 5*sg); h = zeros(size(tg));
      for e = 1:Nel, tt = tg - te(e); m = abs(tt) < 5*sg; h(m) = h(m) + amp(e)*exp(-tt(m).^2/(2*sg^2)).*exp(2i*pi*f0*tt(m)); end
      tau = (X*s + Z*co)/c; L = numel(tg); hc = h(:);
      for r = 1:Nel
        msk = abs(xe(r) - X) <= Z/(2*Fnum); if ~any(msk(:)), continue; end
        tq = tau(msk) + (sqrt((X(msk) - xe(r)).^2 + Z(msk).^2) - rT(r))/c;
        q = (tq - tg(1))/dt; k = floor(q); fr = q - k; ok = (k >= 0) & (k < L - 1); k(~ok) = 0;
        TW(msk) = TW(msk) + ET(r)*((hc(k + 1).*(1 - fr) + hc(k + 2).*fr).*ok);
      end
    end
    Aimg = abs(TW)/max(abs(TW(:)));
    lev(pass) = 20*log10(max(Aimg(Z > zt + 1)));      % strongest response more than 1 mm below the target
    if pass == 1, img{ic} = 20*log10(Aimg + 1e-6); end
  end
  fprintf('  %4.1f  %4.2f   %3d    %5.3f    %6.1f   %6.2f       %6.1f          %6.1f\n', f0, d/lam, th, 1/(1 + abs(sind(th))), thg, Dl, lev(1), lev(2));
end

%% Figures
figure(1);
nm = {'Tx-only', 'Rx-only', 'two-way'};
for k = 1:3
  subplot(3, 1, k);
  plot(xs, 20*log10(prof{1}{k} + 1e-6), 'k--', xs, 20*log10(prof{3}{k} + 1e-6), 'b-');
  axis([-10 10 -60 3]); grid on; ylabel('dB');
  title([nm{k} ' lateral pattern at z = 30 mm: 1 angle (dashed), 11 angles over +-15 deg']);
end
xlabel('x (mm)');
figure(2);
tl = {'5 MHz, d 0.15 mm (0.49\lambda), 15 deg, z 15 mm', '5 MHz, d 0.30 mm (0.97\lambda), 15 deg, z 15 mm', ...
      '7.5 MHz, d 0.30 mm (1.46\lambda), 0 deg, z 30 mm', '7.5 MHz, d 0.30 mm (1.46\lambda), 15 deg, z 30 mm'};
for ic = 1:4
  subplot(2, 2, ic); zt = cases(ic, 4);
  imagesc(-10:0.15:10, (zt - 3):0.1:(zt + 9), img{ic}, [-60 0]); axis image; colormap(gray); colorbar;
  xlabel('x (mm)'); ylabel('z (mm)'); title(tl{ic});
end
