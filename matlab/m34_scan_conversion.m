%% m34_scan_conversion.m  --  sector scan conversion: line density, interpolation domain, scalloping
% Purpose : A phased-array sector (2.5 MHz, 60 % bandwidth, aperture D = 20 mm, +-40 deg,
%           10-90 mm) is simulated on a fine polar grid as complex IQ speckle (k-space filtering:
%           Gaussian axial band, triangular angular transfer of half-width D/lambda cycles/rad)
%           plus point targets. N beam lines are taken from it and scan converted to a
%           Cartesian grid by bilinear interpolation in (r, theta) in one of three domains:
%             IQ (complex, detect afterwards), envelope |z|, or log 20 log10 |z|.
%           Sampling theory: the two-way IQ field occupies |nu| <= D/lambda (cycles/rad), so the
%           IQ Nyquist line spacing is lambda/(2D); detection doubles the band, lambda/(4D).
%           Metrics in a deep speckle region: mean dB on beam lines vs midway between lines
%           (re the fine-grid mean intensity), log-speckle std, and the worst-case peak loss of
%           a point target lying midway between two lines.
% Change  : Nlist (line counts), D, f0, sector angle.
% Expect  : lambda/(2D) = 0.88 deg (91 lines over 80 deg), lambda/(4D) = 0.44 deg (182 lines).
%           Midway between 48 lines (1.7 deg, twice the IQ Nyquist spacing) IQ interpolation darkens
%           the speckle by about 1.9 dB (radial dark bands), envelope interpolation brightens it by
%           about 0.5 dB and log interpolation keeps the mean dB; the log-speckle spread midway is
%           5.5 dB (IQ), 3.3 dB (envelope), 3.6 dB (log). At 192 lines the on-line/midway
%           differences fall to 0.2 dB. A point target midway between lines loses 7.4 dB (48 lines),
%           1.7 dB (96) and 0.4 dB (192) in every domain.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(3);
c = 1540;  f0 = 2.5e6;  lam = c/f0;  B = 0.6*f0;  D = 20e-3;
th_span = 80*pi/180;  r0 = 10e-3;  r1 = 90e-3;
Nr = 512;  Nt = 1024;                                        % fine polar grid
dr = (r1 - r0)/Nr;  dt = th_span/Nt;
rr = r0 + ((0:Nr-1)' + 0.5)*dr;  tt = -th_span/2 + ((0:Nt-1) + 0.5)*dt;
fprintf('lambda = %.3f mm; two-way -6 dB beam width 0.886 lambda/D = %.2f deg\n', lam*1e3, 0.886*lam/D*180/pi);
fprintf('IQ Nyquist line spacing lambda/(2D) = %.2f deg (%d lines over %.0f deg); detected lambda/(4D) = %.2f deg (%d lines)\n', ...
        lam/(2*D)*180/pi, ceil(th_span/(lam/(2*D))), th_span*180/pi, lam/(4*D)*180/pi, ceil(th_span/(lam/(4*D))));

sinc_two_way = @(dth, lam, D) (sin(pi*D*sin(dth)/lam + eps)./(pi*D*sin(dth)/lam + eps)).^2;   % two-way rect amplitude
ur = [0:Nr/2-1, -Nr/2:-1]'/(Nr*dr);  ut = [0:Nt/2-1, -Nt/2:-1]/(Nt*dt);
sf = B/(2*sqrt(2*log(2)));
H = exp(-(c*ur/2).^2/(2*sf^2)) * max(0, 1 - abs(ut)*lam/D);
W = (randn(Nr, Nt) + 1i*randn(Nr, Nt))/sqrt(2);
h = ifft2(H);  apt = sqrt(3000*sum(abs(h(:)).^2)/max(abs(h(:)).^2));
Zs = ifft2(fft2(W) .* H);                                   % speckle only, for the reference level
pt_r = [25e-3 40e-3 55e-3];  pt_t = [0 5 -10]*pi/180;
for q = 1:3
    [~, ir] = min(abs(rr - pt_r(q)));  [~, it] = min(abs(tt - pt_t(q)));  W(ir, it) = W(ir, it) + apt;
end
Zf = ifft2(fft2(W) .* H);
spk = rr > 60e-3 & rr < 85e-3;  tsp = abs(tt) < 30*pi/180;
I0 = mean(mean(abs(Zs(spk, tsp)).^2));                       % fine-grid mean speckle intensity

%% Cartesian output grid
Nx = 480;  Nz = 360;  xg = linspace(-60e-3, 60e-3, Nx);  zg = linspace(0, 90e-3, Nz);
[XG, ZG] = meshgrid(xg, zg);  RG = hypot(XG, ZG);  TG = atan2(XG, ZG);
inside = RG >= r0 & RG <= rr(end) & abs(TG) <= th_span/2 - 1e-9;
fr = (RG - rr(1))/dr + 1;                                    % fractional sample index (1-based)
regS = inside & RG > 60e-3 & RG < 85e-3 & abs(TG) < 30*pi/180;

Nlist = [48 96 192];  doms = {'IQ', 'env', 'log'};
fprintf('\n lines  spacing  domain | mean dB on-line  midway | log std on-line  midway | worst point loss (dB)\n');
for N = Nlist
    tl = linspace(-th_span/2 + dt/2, th_span/2 - dt/2, N);   % beam line angles
    il = round((tl + th_span/2)/dt + 0.5);  il = min(Nt, max(1, il));
    Zl = Zf(:, il);  tl = tt(il);
    ft = interp1(tl, 1:N, TG, 'linear', NaN);                % fractional line index
    dlt = abs(ft - round(ft));  onl = regS & dlt < 0.1;  mid = regS & dlt > 0.4;
    for q = 1:3
        switch doms{q}
            case 'IQ',  V = Zl;
            case 'env', V = abs(Zl);
            case 'log', V = 20*log10(abs(Zl));
        end
        i0 = floor(fr);  a = fr - i0;  j0 = floor(ft);  b = ft - j0;
        ok = inside & ~isnan(ft) & i0 >= 1 & i0 < Nr & j0 >= 1 & j0 < N;
        j0(~ok) = 1;  i0(~ok) = 1;  a(~ok) = 0;  b(~ok) = 0;
        j1 = min(N, j0 + 1);
        idx = @(i, j) i + (j - 1)*Nr;
        Vi = (1 - a).*(1 - b).*V(idx(i0, j0)) + a.*(1 - b).*V(idx(i0 + 1, j0)) + (1 - a).*b.*V(idx(i0, j1)) + a.*b.*V(idx(i0 + 1, j1));
        if q == 1, Dimg = 10*log10(abs(Vi).^2/I0); elseif q == 2, Dimg = 20*log10(Vi/sqrt(I0)); else, Dimg = Vi - 10*log10(I0); end
        Dimg(~ok) = -60;
        % point target midway between two lines: amplitude of a point at angle tm sampled by lines at tm -/+ half spacing
        pl = 20*log10(sinc_two_way(0.5*(tl(2) - tl(1)), lam, D));
        fprintf(' %4d   %5.2f deg  %-4s  |     %6.2f       %6.2f  |      %5.2f      %5.2f  |   %5.1f\n', N, (tl(2) - tl(1))*180/pi, doms{q}, ...
                mean(Dimg(onl)), mean(Dimg(mid)), std(Dimg(onl)), std(Dimg(mid)), pl);
        if N == 48 && q == 1, img48 = Dimg; end
        if N == 192 && q == 1, img192 = Dimg; end
    end
end
fprintf('(fully developed speckle: mean dB of a pixel = -2.51 dB, std 5.57 dB re 10 log10 <I>)\n');

figure('Position', [60 60 1100 420]);
subplot(1, 2, 1); imagesc(xg*1e3, zg*1e3, img48, [-40 15]); axis image; colormap(gray(256));
title('48 lines, IQ-domain bilinear'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(1, 2, 2); imagesc(xg*1e3, zg*1e3, img192, [-40 15]); axis image;
title('192 lines, IQ-domain bilinear'); xlabel('x (mm)');

