%% m21_spatial_impulse_response.m
% Spatial impulse response (SIR) of a baffled rectangular element, exact and Field II-style.
% Everything here is computed by this script; Field II itself is not called (the subdivision
%   only reimplements its far-field sub-element idea).
%   h(r,t) = (c/2pi) * Theta(t), Theta = angle of the circle of radius sqrt((ct)^2 - z^2)
%   (centred on the projection of r) that lies inside the aperture.  p = rho * dv/dt * h.
% What it shows
%   1) exact SIR on axis, inside the edge and beyond the edge; check against the closed form
%   2) integral of h equals the surface integral of 1/(2 pi R)
%   3) Field II-style subdivision (point-like sub-elements) converging to the exact SIR
%   4) on-axis CW field of a 6 mm square: subdivided sum versus the exact SIR
% What to change: element size (w, H), field points, f0, sub-element counts.
% Expected output: plateau h/c = 1.000 on axis, off-axis duration 1.129 us and peak 0.204c,
%   last on-axis CW maximum of the 6 mm square near 24 mm, errors ~22, 4.4, 0.6, 0.13 dB.
% Units: mm, us, MHz; c in mm/us.  Base MATLAB / Octave.

clear; close all; rng(1);
c   = 1.54;            % mm/us (1540 m/s)
w   = 4;  H = 6;       % element width (x) and height (y), mm
a   = w/2; b = H/2;
fs  = 200;             % MHz, SIR sampling
nAng = 1440;           % angles on the arrival circle
phi = ((0:nAng-1) + 0.5) * 2*pi/nAng;

% exact SIR by testing points on the arrival circle (vectorized over angles)
sirRect = @(x0, y0, z0, t) c * mean( ...
    (abs(x0 + sqrt(max((c*t(:)).^2 - z0^2, 0)) * cos(phi)) < a) & ...
    (abs(y0 + sqrt(max((c*t(:)).^2 - z0^2, 0)) * sin(phi)) < b) & ...
    repmat((c*t(:)).^2 > z0^2, 1, nAng), 2).';     % returns a row vector

%% 1) on axis, inside the edge, beyond the edge (z = 15 mm)
z0  = 15;
pts = [0 0; 1.5 0; 6 0];
names = {'on axis', 'x = 1.5 mm (inside)', 'x = 6 mm (beyond edge)'};
f0 = 3; B = 0.6;                                   % pulse: 3 MHz, 60 % -6 dB bandwidth
sf = B*f0/(2*sqrt(2*log(2))); st = 1/(2*pi*sf);
dv = @(t) exp(-t.^2/(2*st^2)) .* (2*pi*f0*cos(2*pi*f0*t) - t/st^2 .* sin(2*pi*f0*t));
tk = (-ceil(4*st*fs):ceil(4*st*fs))/fs;            % pulse kernel time axis
figure(1); clf;
for k = 1:3
    x0 = pts(k,1); y0 = pts(k,2);
    Rmin = sqrt(max(abs(x0)-a,0)^2 + max(abs(y0)-b,0)^2 + z0^2);
    Rmax = sqrt((abs(x0)+a)^2 + (abs(y0)+b)^2 + z0^2);
    t = (floor(Rmin/c*fs)-2 : ceil(Rmax/c*fs)+2)/fs;
    h = sirRect(x0, y0, z0, t);
    corners = sort(sqrt((x0-[-a a -a a]).^2 + (y0-[-b -b b b]).^2 + z0^2)/c);
    fprintf('%-24s t0 = %.3f us, t1 = %.3f us, duration = %.3f us, peak h/c = %.3f\n', ...
        names{k}, Rmin/c, Rmax/c, (Rmax-Rmin)/c, max(h)/c);
    fprintf('%-24s corner arrivals: %s us\n', '', sprintf('%.3f ', corners));
    p = conv(h, dv(tk)) / fs;                       % rho = 1
    tp = t(1) + tk(1) + (0:numel(p)-1)/fs;
    pff = w*H/(2*pi*sqrt(x0^2+y0^2+z0^2)) * max(abs(dv(tk)));
    fprintf('%-24s pressure peak / point-source estimate = %.3f\n', '', max(abs(p))/pff);
    subplot(2,1,1); plot(t - z0/c, h/c, 'LineWidth', 1.5); hold on;
    subplot(2,1,2); plot(tp - z0/c, p/pff, 'LineWidth', 1.2); hold on;
    if k == 1
        % closed form on axis: Theta = 2pi - 4acos(a/s)[s>a] - 4acos(b/s)[s>b]
        s = sqrt(max((c*t).^2 - z0^2, 0));
        th = 2*pi*ones(size(s));
        th(s > a) = th(s > a) - 4*acos(a./s(s > a));
        th(s > b) = th(s > b) - 4*acos(b./s(s > b));
        th(s >= sqrt(a^2+b^2) | (c*t).^2 <= z0^2) = 0;
        fprintf('%-24s max |numeric - closed form| / c = %.1e\n', '', max(abs(h - c*th/(2*pi)))/c);
    end
end
subplot(2,1,1); ylim([0 1.15]); xlabel('t - z/c (us)'); ylabel('h / c'); title('Spatial impulse response, 4 x 6 mm, z = 15 mm'); legend(names); grid on;
subplot(2,1,2); xlabel('t - z/c (us)'); ylabel('p / point-source peak'); title('Pressure p = \rho v'' * h, 3 MHz 60 % pulse'); grid on;

%% 2) integral check: int h dt = int dS / (2 pi R)
x0 = 5; y0 = 1; z0 = 8;
Rmin = sqrt(max(abs(x0)-a,0)^2 + max(abs(y0)-b,0)^2 + z0^2); Rmax = sqrt((abs(x0)+a)^2 + (abs(y0)+b)^2 + z0^2);
t = (floor(Rmin/c*800):ceil(Rmax/c*800))/800;
I1 = sum(sirRect(x0, y0, z0, t)) / 800;
ng = 400; [XS, YS] = meshgrid(-a + ((1:ng)-0.5)*w/ng, -b + ((1:ng)-0.5)*H/ng);
I2 = sum(sum((w/ng)*(H/ng) ./ (2*pi*sqrt((x0-XS).^2 + (y0-YS).^2 + z0^2))));
fprintf('\nIntegral check at (5,1,8) mm: int h dt = %.5f mm, surface integral = %.5f mm\n', I1, I2);

%% 3) Field II-style subdivision of the SIR at the off-axis point (6, 0, 15)
x0 = 6; y0 = 0; z0 = 15; fsb = 200;
Rmin = sqrt((x0-a)^2 + z0^2); Rmax = sqrt((x0+a)^2 + b^2 + z0^2);
t = (floor(Rmin/c*fsb)-4 : ceil(Rmax/c*fsb)+4)/fsb;
hx = sirRect(x0, y0, z0, t);
pex = conv(hx, dv(tk))/fsb;
fprintf('\nSubdivision at (6, 0, 15) mm, pressure-pulse error vs exact:\n');
for n = [2 4 8 16 32]
    [XS, YS] = meshgrid(-a + ((1:n)-0.5)*w/n, -b + ((1:n)-0.5)*H/n);
    R = sqrt((x0-XS(:)).^2 + (y0-YS(:)).^2 + z0^2); amp = (w/n)*(H/n) ./ (2*pi*R);
    u = (R/c - t(1))*fsb; i0 = floor(u); fr = u - i0;                       % linear split of each delta
    hs = accumarray([i0+1; i0+2], [amp.*(1-fr); amp.*fr], [numel(t)+2, 1])' * fsb;
    hs = hs(1:numel(t));
    ps = conv(hs, dv(tk))/fsb;
    fprintf('  %2d x %2d sub-elements (%.2f x %.2f mm): rms error %.1f %% of peak\n', n, n, w/n, H/n, 100*sqrt(mean((ps-pex).^2))/max(abs(pex)));
end

%% 4) on-axis CW field of a 6 mm square at 3 MHz: exact SIR vs subdivision
ws = 6; as = 3; f = 3; k = 2*pi*f/c; lam = c/f;
zs = 2:0.25:60;
Pex = zeros(size(zs));
for iz = 1:numel(zs)
    z = zs(iz); t0 = z/c; t1 = sqrt(z^2 + 2*as^2)/c; nt = max(400, ceil((t1-t0)*2000));
    tt = t0 + ((1:nt)-0.5)*(t1-t0)/nt; s = sqrt((c*tt).^2 - z^2);
    th = 2*pi - 8*acos(min(as./s, 1));                                     % square: a = b
    th(s >= sqrt(2)*as) = 0;
    Pex(iz) = abs(sum(c*th/(2*pi) .* exp(-1i*2*pi*f*tt)) * (t1-t0)/nt);
end
mx = false(size(Pex)); WN = 24;
for i = WN+1:numel(Pex)-WN
    mx(i) = Pex(i) == max(Pex(i-WN:i+WN));
end
iz = find(mx, 1, 'last');
fprintf('\n6 mm square, 3 MHz: last on-axis CW maximum at z = %.2f mm (a^2/lambda = %.2f mm, ratio %.2f)\n', zs(iz), as^2/lam, zs(iz)/(as^2/lam));
figure(2); clf; plot(zs, 20*log10(Pex/max(Pex)), 'k', 'LineWidth', 1.8); hold on;
for n = [1 8 16 32]
    [XS, YS] = meshgrid(-as + ((1:n)-0.5)*ws/n, -as + ((1:n)-0.5)*ws/n); dA = (ws/n)^2;
    R = sqrt(XS(:).^2 + YS(:).^2 + zs.^2);                                  % (n^2) x numel(zs)
    Pap = abs(sum(dA ./ (2*pi*R) .* exp(-1i*k*R), 1));
    e = max(abs(20*log10(Pap(zs >= 5) ./ Pex(zs >= 5))));
    fprintf('  %2d x %2d sub-elements (%.3f mm = %.2f lambda): max on-axis error (z >= 5 mm) = %.2f dB\n', n, n, ws/n, ws/n/lam, e);
    plot(zs, 20*log10(Pap/max(Pex)), '--', 'LineWidth', 1.1);
end
xlabel('z (mm)'); ylabel('|p| on axis (dB re max)'); ylim([-30 3]); grid on;
title('6 mm square, 3 MHz CW: exact SIR (black) vs subdivided sums');
legend('exact', '1x1', '8x8', '16x16', '32x32', 'Location', 'southeast');
