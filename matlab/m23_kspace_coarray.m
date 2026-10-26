%% m23_kspace_coarray.m  --  the coarray: two-way PSF = Fourier transform of (Tx aperture * Rx aperture)
% Purpose : narrowband focal-plane (paraxial) model at depth z. A transmit "aperture" (real
%           elements, or virtual sources xi_n = z sin(theta_n) for plane waves) and a receive
%           aperture give one-way factors T(x), R(x). The two-way PSF T(x)R(x) is the Fourier
%           transform of the coarray: every Tx/Rx pair (xi_t, xi_r) adds a_t a_r at the
%           spatial frequency kx = (xi_t + xi_r)/(lambda z). The script checks that identity,
%           then compares designs with the same 64-element receive aperture: full, Hann Tx,
%           one element (or one plane wave), sparse periodic, vernier (Tx every 3rd, Rx every
%           2nd), random thinning, and plane-wave compounding (3, 11, 31 angles).
% Change  : z, the cases, the random thinning fraction, the plane-wave angle sets.
% Expect  : identity error ~1e-13. Two-way PSL: full -26.5 dB; single element / one plane
%           wave -13.3 dB (receive only); sparse every-2nd: grating lobe at
%           x = lambda z/(2d) = 15.4 mm at ~0 dB; vernier 3/2: no grating lobe, PSL ~ -27 dB;
%           plane waves 1 deg step: replicas at lambda/d(sin) = 17.6 mm, beyond the window.
%           One plane wave or one element: coarray = receive aperture only. 31 angles at 1 deg:
%           the replica at 17.9 mm survives at the receive sidelobe level there (-38 dB).
%           Random thinning: mean one-way sidelobe power (N-K)/(K(N-1)) within 0.7 dB.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

c = 1.540; f0 = 5; lam = c/f0;        % mm/us, MHz, mm
N = 64; d = 0.30; z = 30;              % elements, pitch (mm), depth (mm)
xe = ((0:N-1) - (N-1)/2) * d;          % element positions (mm)
x  = linspace(-20, 20, 8001).';        % lateral position at depth z (mm)
FT = @(pos, w) exp(-1i*2*pi*x*pos(:).'/(lam*z)) * w(:);   % one-way factor, paraxial
fprintf('N = %d, pitch %.2f mm (%.2f lambda), z = %g mm, aperture %.1f mm (F# %.2f)\n', ...
    N, d, d/lam, z, N*d, z/(N*d));

%% (1) Identity: FT of the coarray equals the product of the one-way factors
at = 0.5 - 0.5*cos(2*pi*((0:N-1) + 0.5)/N);   ar = ones(1, N);       % Hann Tx, rect Rx
co = conv(at, ar);                                                    % coarray, lags in pitch units
pos2 = ((0:2*N-2) - (N-1)) * d;                                       % xi_t + xi_r (mm)
lhs = FT(pos2, co);  rhs = FT(xe, at) .* FT(xe, ar);
fprintf('(1) max |FT(coarray) - T R| / max|T R| = %.1e\n\n', max(abs(lhs - rhs))/max(abs(rhs)));

%% (2) Designs with the same receive aperture
th3  = (-15:15:15)*pi/180;  th11 = (-15:3:15)*pi/180;  th31 = (-15:1:15)*pi/180;
evr2 = mod(0:N-1, 2) == 0;  evr3 = mod(0:N-1, 3) == 0;
rnd  = @(K) full(sparse(1, randperm(N, K), 1, 1, N));
cases = { ...
  'full Tx = Rx (focused)',          xe,            ones(1,N),       xe, ones(1,N);
  'Hann Tx, rect Rx',                xe,            at,              xe, ones(1,N);
  'one element / one plane wave',    0,             1,               xe, ones(1,N);
  'sparse: Tx = Rx every 2nd',       xe(evr2),      ones(1,sum(evr2)), xe(evr2), ones(1,sum(evr2));
  'vernier: Tx every 3rd, Rx 2nd',   xe(evr3),      ones(1,sum(evr3)), xe(evr2), ones(1,sum(evr2));
  'random 50 %, Tx ~= Rx',           xe,            rnd(N/2),        xe, rnd(N/2);
  'plane waves 3 x 15 deg',          z*sin(th3),    ones(1,3),       xe, ones(1,N);
  'plane waves 11 x 3 deg',          z*sin(th11),   ones(1,11),      xe, ones(1,N);
  'plane waves 31 x 1 deg',          z*sin(th31),   ones(1,31),      xe, ones(1,N)};
nc = size(cases, 1);  P = zeros(numel(x), nc);
fprintf('(2) %-32s %4s %4s %8s %7s %7s %8s %9s %9s\n', 'design', 'nTx', 'nRx', 'coarray', 'w-6dB', 'PSL', 'at x', 'GL/rep x', 'level');
fprintf('    %-32s %4s %4s %8s %7s %7s %8s %9s %9s\n', '', '', '', '(mm)', '(mm)', '(dB)', '(mm)', '(mm)', '(dB)');
xg = [NaN NaN NaN lam*z/(2*d) NaN NaN lam/mean(diff(sin(th3))) lam/mean(diff(sin(th11))) lam/mean(diff(sin(th31)))];
mid = (numel(x)+1)/2;
for k = 1:nc
    [pt, wt, pr, wr] = cases{k, 2:5};
    wt = wt(:).'; wr = wr(:).'; keepT = wt ~= 0; keepR = wr ~= 0;
    pt = pt(keepT); wt = wt(keepT); pr = pr(keepR); wr = wr(keepR);
    p2 = abs(FT(pt, wt) .* FT(pr, wr));  p2 = p2/p2(mid);  P(:,k) = p2;
    span = (max(pt) + max(pr)) - (min(pt) + min(pr));           % coarray extent (mm)
    iN = mid - 1 + find(diff(p2(mid:end)) > 0, 1);              % first null right of center
    w6 = 2*x(mid - 1 + find(p2(mid:end) < 0.5, 1));
    sl = p2;  sl(2*mid - iN:iN) = 0;  [pk, ipk] = max(sl);
    lev = NaN; if ~isnan(xg(k)), lev = 20*log10(max(p2(abs(x - xg(k)) <= 0.5))); end
    fprintf('    %-32s %4d %4d %8.1f %7.2f %7.1f %8.1f %9.2f %9.1f\n', cases{k,1}, numel(pt), numel(pr), span, w6, 20*log10(pk), abs(x(ipk)), xg(k), lev);
end
fprintf('    predicted: sparse GL at lambda z/(2d) = %.1f mm; plane-wave replicas lambda/d(sin): %.2f, %.2f, %.2f mm\n', ...
    lam*z/(2*d), lam/mean(diff(sin(th3))), lam/mean(diff(sin(th11))), lam/mean(diff(sin(th31))));

%% (3) Random thinning: mean sidelobe floor vs number of active elements K
fprintf('\n(3) random thinning, one-way, 300 trials each\n%6s %12s %18s\n', 'K', 'mean SL (dB)', '(N-K)/(K(N-1))');
s  = linspace(-1, 1, 4001).';  far = abs(s) > 4/N;      % sin(theta), outside the main lobe
E  = exp(1i*pi*s*(0:N-1));                               % steering matrix, d = lambda/2
for K = [8 16 32 48]
    acc = 0;
    for tr = 1:300
        w = rnd(K);  AF = abs(E * w(:)).^2 / K^2;
        acc = acc + mean(AF(far));
    end
    fprintf('%6d %12.2f %18.2f\n', K, 10*log10(acc/300), 10*log10((N-K)/(K*(N-1))));
end

%% Figure: coarrays and two-way PSFs
show = [1 3 4 5 8 9];
figure('Position', [40 40 1200 700]);
for j = 1:numel(show)
    k = show(j); [pt, wt, pr, wr] = cases{k, 2:5};
    % coarray binned on a grid of one pitch in xi_t + xi_r (virtual sources fall between)
    g  = -40:d:40;  cg = zeros(size(g));
    for a = 1:numel(pt)
        if wt(a) == 0, continue; end
        idx = round((pt(a) + pr - g(1))/d) + 1;
        cg  = cg + accumarray(idx(:), wt(a)*wr(:), [numel(g) 1]).';
    end
    subplot(3, 4, 2*j - 1); stem(g(cg ~= 0), cg(cg ~= 0), 'Marker', 'none'); xlim([-40 40]); grid on;
    title(cases{k,1}); xlabel('\xi_t + \xi_r (mm)'); ylabel('coarray weight');
    subplot(3, 4, 2*j); plot(x, 20*log10(P(:,k) + 1e-6), 'LineWidth', 1.1); ylim([-60 3]); grid on;
    xlabel('x at z = 30 mm (mm)'); ylabel('two-way (dB)');
end
