%% m31_saturation_frc.m
% Purpose : two ways to ask "is my ULM image finished, and how sharp is it?"
%   Part 1 - saturation. Time to fill a vessel of diameter d with pixels of
%           size l (Hingot et al., Sci. Rep. 2019, Eq. 4 and 6):
%               T(d) = max(d, l) / (eta * C * Q(d) * l),  Q = pi d^2 v / 4,
%           i.e. each bubble paints one pixel column of the max(d,l)/l columns
%           across the vessel. Saturation of that vessel: S(t) = 1 - exp(-t/T).
%           Checked by a Poisson Monte Carlo with uniform and with flux-weighted
%           (parabolic) lateral bubble positions.
%   Part 2 - Fourier ring correlation (FRC). Localizations along a random
%           vessel network are split into two halves, rendered, and their
%           ring-wise spectral correlation is read against three thresholds:
%           fixed 1/7, half-bit (van Heel and Schatz 2005) and 2-sigma.
%           Split by track (odd/even tracks, as in Hingot et al. 2021) versus
%           split by localization (odd/even points of every track).
% Model   : v = kv * d (kv = 0.2 mm/s per um: 8 um -> 1.6 mm/s, 50 um -> 10 mm/s),
%           C = 200 bubbles per uL (2e5 per mL), detection efficiency eta.
%           FRC: 2.56 x 2.56 mm field, 5 um render pixels (512 x 512); 40
%           vessels of 10-80 um, tracks allocated in proportion to flux (d^3),
%           points every 10 um along a fixed streamline, Gaussian errors sigma.
% Change  : C, kv, eta, lpix; nTrackList, sigList.
% Expect  : T = 255 s for 5 um capillaries at any pixel size, 15.9 s for a
%           20 um vessel with 5 um pixels, 0.6 s for 100 um (eta = 1). The
%           Monte Carlo matches 1 - exp(-t/T) for uniform lateral positions
%           (0.634 vs 0.632 at t = T); a parabolic profile fills the wall
%           pixels later (0.862 instead of 0.900 at 2.3 T). FRC: track-split
%           resolution improves from 64 to 36 um (sigma 5 um, 1/7) as tracks
%           grow from 300 to 10,000: it measures sampling as well as precision.
%           Splitting by localization puts the same streamline in both halves
%           and reports 28-31 um from the start: optimistic when undersampled.
%           Half-bit is the most conservative threshold at low counts.
% Runtime : about 20 s in Octave.
rng(1);

%% Part 1: saturation times
C = 200;            % bubbles per microlitre (2e5 per mL)
kv = 0.2;           % mm/s per um of diameter
eta = 1;            % fraction of passing bubbles that end up in kept tracks
dList = [5 10 20 50 100];  lList = [5 10 20];
Tsat = @(d, l) max(d, l)./(eta*C*(pi*(d*1e-6).^2/4.*(kv*d*1e-3))*1e9 .* l);   % C (/uL) x Q (uL/s; 1 m^3 = 1e9 uL)
fprintf('Saturation time constant T (s) = max(d,l)/(eta C Q l), C = %d /uL, v = %.1f mm/s per um, eta = %.2f\n', C, kv, eta);
fprintf('   d (um)   v (mm/s)   Q (nL/s)   %s\n', sprintf('   l = %2d um', lList));
for d = dList
    Q = pi*(d*1e-6)^2/4*(kv*d*1e-3)*1e12;  % nL/s
    fprintf('   %5d    %6.1f    %8.4f   %s\n', d, kv*d, Q, sprintf('%10.1f  ', Tsat(d, lList)));
end
fprintf('   time to 90 %% = T ln 10 = 2.30 T; 3D (voxel l^3) multiplies T by max(d,l)/l (Hingot Eq. 4'')\n');

% Monte Carlo check for d = 20 um, l = 5 um (4 pixel columns)
d = 20; l = 5; T = Tsat(d, l); ncol = d/l; tgrid = linspace(0, 4*T, 41); nRep = 2000;
rate = eta*C*pi*(d*1e-6)^2/4*(kv*d*1e-3)*1e9;           % bubbles per second through the vessel
Su = zeros(size(tgrid)); Sp = Su;
for r = 1:nRep
    tb = cumsum(-log(rand(1, ceil(rate*tgrid(end)*3) + 10)))/rate; tb = tb(tb <= tgrid(end));
    cu = ceil(rand(size(tb))*ncol);                        % uniform across the lumen
    up = zeros(size(tb)); for k = 1:numel(tb), u = 2*rand - 1; while rand > 1 - u^2, u = 2*rand - 1; end; up(k) = u; end
    cp = min(ncol, floor((up + 1)/2*ncol) + 1);            % flux-weighted (parabolic profile)
    for it = 1:numel(tgrid)
        Su(it) = Su(it) + numel(unique(cu(tb <= tgrid(it))))/ncol;
        Sp(it) = Sp(it) + numel(unique(cp(tb <= tgrid(it))))/ncol;
    end
end
Su = Su/nRep; Sp = Sp/nRep;
i1 = find(tgrid >= T, 1); i23 = find(tgrid >= 2.3*T, 1);
fprintf('\nMonte Carlo, d = %d um, l = %d um (T = %.1f s): S(T) = %.3f uniform, %.3f flux-weighted (model %.3f)\n', d, l, T, Su(i1), Sp(i1), 1 - exp(-1));
fprintf('                                     S(2.3T) = %.3f uniform, %.3f flux-weighted (model %.3f)\n', Su(i23), Sp(i23), 0.9);

%% Part 2: FRC on a simulated vessel network
Np = 512; pr = 5; Lf = Np*pr;                              % render pixels (um), field (um)
nV = 40; dv = exp(log(10) + rand(1,nV)*log(8));            % diameters 10-80 um, log-uniform
VX = cell(1,nV); VZ = VX; VN = VX;
for k = 1:nV                                               % smooth random curves across the field
    th = rand*pi; c = rand(1,2)*Lf; s = linspace(-Lf, Lf, 4000)';
    a = 150*randn; f = 1/(400 + 800*rand);
    x = c(1) + s*cos(th) - a*sin(2*pi*f*s)*sin(th); z = c(2) + s*sin(th) + a*sin(2*pi*f*s)*cos(th);
    in = x > 0 & x < Lf & z > 0 & z < Lf; x = x(in); z = z(in);
    if numel(x) < 50, x = c(1) + [0;1]; z = c(2) + [0;1]; end
    cum = [0; cumsum(sqrt(diff(x).^2 + diff(z).^2))]; [cum, iu] = unique(cum); x = x(iu); z = z(iu);
    su = (0:2:cum(end))'; VX{k} = interp1(cum, x, su); VZ{k} = interp1(cum, z, su);
    gx = gradient(VX{k}); gz = gradient(VZ{k}); gl = sqrt(gx.^2 + gz.^2) + eps; VN{k} = [-gz./gl gx./gl];
end
flux = dv.^3; cwf = cumsum(flux)/sum(flux);
q = (0:Np-1); q(q >= Np/2) = q(q >= Np/2) - Np; [QX, QZ] = meshgrid(q, q); QR = sqrt(QX.^2 + QZ.^2);
ring = round(QR); nR = Np/2; nring = accumarray(ring(:) + 1, 1); nring = nring(1:nR);
thrHalf = (0.2071 + 1.9102./sqrt(nring))./(1.2071 + 0.9102./sqrt(nring));
thr2s = 2./sqrt(nring/2);
qf = (0:nR-1)'/(Np*pr);                                    % cycles per um
nTrackList = [300 1000 3000 10000]; sigList = [5 10 20];
fprintf('\nFRC resolution (um) = 1/q at the first crossing; render pixel %d um, Nyquist limit %d um\n', pr, 2*pr);
fprintf('   sigma  tracks   locs     split     1/7     half-bit   2-sigma\n');
keepCurves = {};
for sg = sigList
    for nT = nTrackList
        P1 = cell(1, nT); nloc = 0;                                                     % [iz ix track-parity point-parity]
        for t = 1:nT
            k = find(rand <= cwf, 1); n = numel(VX{k}); len = round((200 + 400*rand)/2);   % 0.2-0.6 mm of track, points every 2 um of curve
            if n < len + 2, len = n - 2; end
            j0 = 1 + floor(rand*(n - len)); jj = (j0:5:j0+len)';                        % a point every 10 um
            u = 2*rand - 1; while rand > 1 - u^2, u = 2*rand - 1; end                   % fixed streamline (flux-weighted)
            px = VX{k}(jj) + u*dv(k)/2*VN{k}(jj,1) + sg*randn(numel(jj),1);
            pz = VZ{k}(jj) + u*dv(k)/2*VN{k}(jj,2) + sg*randn(numel(jj),1);
            ix = floor(px/pr) + 1; iz = floor(pz/pr) + 1; ok = ix >= 1 & ix <= Np & iz >= 1 & iz <= Np;
            P1{t} = [iz(ok) ix(ok) mod(t,2)*ones(sum(ok),1) mod((1:sum(ok))',2)];
        end
        PP = cat(1, P1{:}); nloc = size(PP,1);
        im1 = accumarray(PP(PP(:,3)==1,1:2), 1, [Np Np]); im2 = accumarray(PP(PP(:,3)==0,1:2), 1, [Np Np]);     % split by track
        im1l = accumarray(PP(PP(:,4)==1,1:2), 1, [Np Np]); im2l = accumarray(PP(PP(:,4)==0,1:2), 1, [Np Np]);   % split by localization
        for sp = 1:2
            if sp == 1, F1 = fft2(im1); F2 = fft2(im2); lab = 'tracks'; else, F1 = fft2(im1l); F2 = fft2(im2l); lab = 'locs  '; end
            num = accumarray(ring(:) + 1, real(F1(:).*conj(F2(:)))); d1 = accumarray(ring(:) + 1, abs(F1(:)).^2); d2 = accumarray(ring(:) + 1, abs(F2(:)).^2);
            frc = num(1:nR)./sqrt(d1(1:nR).*d2(1:nR)); frc = conv(frc, ones(3,1)/3, 'same'); frc(1) = 1;
            r = zeros(1,3); th = {1/7*ones(nR,1), thrHalf, thr2s};
            for k = 1:3                                    % first crossing above 1/(512 um), i.e. from ring 6
                ic = find(frc(6:end) < th{k}(6:end), 1) + 5;
                if isempty(ic), r(k) = 2*pr; elseif ic == 6, r(k) = Inf; else, r(k) = 1/qf(ic); end
            end
            fprintf('   %4d   %6d  %7d    %s   %6.1f    %6.1f    %6.1f\n', sg, nT, nloc, lab, r);
            if sg == 10 && (nT == 1000 || nT == 10000), keepCurves{end+1} = {frc, sprintf('%d tracks, %s', nT, strtrim(lab))}; end %#ok<SAGROW>
        end
    end
end

%% Figures
figure('Position', [50 50 1100 380]);
subplot(1,3,1); semilogx(dList, Tsat(dList, 5), 'o-', dList, Tsat(dList, 10), 's-', dList, Tsat(dList, 20), 'd-'); grid on;
set(gca, 'YScale', 'log'); xlabel('vessel diameter d (\mum)'); ylabel('T (s)'); legend('l = 5 \mum', 'l = 10 \mum', 'l = 20 \mum'); title('saturation time constant');
subplot(1,3,2); plot(tgrid/T, Su, 'o', tgrid/T, Sp, 's', tgrid/T, 1 - exp(-tgrid/T), 'k-'); grid on;
xlabel('t / T'); ylabel('fraction of pixel columns visited'); legend('MC uniform', 'MC flux-weighted', '1 - e^{-t/T}', 'Location', 'southeast'); title('d = 20 \mum, l = 5 \mum');
subplot(1,3,3); hold on; for k = 1:numel(keepCurves), plot(qf*1000, keepCurves{k}{1}); end
plot(qf*1000, thrHalf, 'k--', qf*1000, thr2s, 'k:', qf*1000, ones(nR,1)/7, 'k-.'); grid on; ylim([-0.1 1]);
legend([cellfun(@(c) c{2}, keepCurves, 'UniformOutput', false), {'half-bit', '2\sigma', '1/7'}]);
xlabel('spatial frequency (cycles/mm)'); ylabel('FRC'); title('FRC, \sigma = 10 \mum');
