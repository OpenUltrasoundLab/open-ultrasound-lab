%% m12_nuisance_fisher.m
% Purpose : what unknown (nuisance) parameters cost, from the Fisher matrix.
%   The bound on parameter 1 with nuisances is [inv(I)](1,1) >= 1/I(1,1).
%   We print the std inflation  sqrt([inv(I)](1,1) * I(1,1))  (1 = free).
%   Noise: circular complex Gaussian, E|w|^2 = 1 per sample, so
%   I = 2 Re(J' * J) with J the Jacobian of the noise-free signal.
%   Part 1: one bubble, Gaussian PSF (FWHM = 1): unknown complex amplitude,
%           width, complex background; centred vs near the window edge.
%   Part 2: a second bubble at separation d (FWHM units), relative phase psi,
%           its position and amplitude unknown. Inflation vs d.
%   Part 3: axial localisation of a carrier-modulated echo: known vs unknown
%           phase. With unknown phase the carrier term cancels exactly.
%   Part 4: sound speed from receive delays of one point target
%           (0 deg plane-wave transmit, 128 elements, pitch 0.3 mm,
%           1 ns independent delay jitter per element).
% Change  : window half-width, offsets, phases, the array, sigma_tau.
% Expect  : Part 1: 1.000 for a centred bubble with every nuisance unknown;
%           ~4.2 when the bubble sits at the edge of a +-1.5 FWHM window.
%           Part 2: in phase 3.96 / 1.39 / 1.00 at d = 0.5 / 1 / 2 FWHM,
%           quadrature 1.81 / 1.11 / 1.00; same values for a weaker neighbour.
%           Part 3: unknown-phase Fisher = envelope-only Fisher (ratio 1.000).
%           Part 4 at z = 30 mm: sigma_c = 0.0034 (z known), 0.065 (z
%           unknown), 2.83 m/s (z and t0 unknown).
rng(1);
fw = 1; w = fw/(2*sqrt(2*log(2))); pix = 1/16;
gfun = @(x, x0) exp(-(x - x0).^2/(2*w^2));
first = @(v) v(1);
inv11 = @(M) first(M \ [1; zeros(size(M, 1) - 1, 1)]);     % [inv(M)](1,1)
infl = @(J) sqrt(inv11(2*real(J'*J)) * 2*real(J(:,1)'*J(:,1)));

%% Part 1: single bubble
fprintf('Part 1: one bubble, std inflation of x1 (1.000 = nuisance is free)\n');
cases = {'centred, window +-3 FWHM', 3, 0; 'offset 1.0, window +-1.5', 1.5, 1.0; 'offset 1.5, window +-1.5 (edge)', 1.5, 1.5};
for k = 1:size(cases, 1)
    win = cases{k, 2}; off = cases{k, 3};
    x = (-round(win/pix):round(win/pix))'*pix; g = gfun(x, off);
    dx = g.*(x - off)/w^2; dw = g.*(x - off).^2/w^3; one = ones(size(x));
    J1 = [dx, g, 1j*g];                       % x1, Re a, Im a
    J2 = [J1, dw];                            % + width
    J3 = [J2, one, 1j*one];                   % + complex background
    fprintf('  %-34s amp: %.3f   amp+width: %.3f   amp+width+bg: %.3f\n', cases{k, 1}, infl(J1), infl(J2), infl(J3));
end

%% Part 2: a neighbour bubble
dList = [0.25 0.5 0.75 1 1.5 2 3];
psiList = [0 90 180]; relList = [1 0.3];
x = (-4/pix:4/pix)'*pix; ratio2 = zeros(numel(psiList), numel(dList));
fprintf('\nPart 2: second bubble at separation d (x2, a2 unknown; a1 unknown)\n');
fprintf('  psi(deg) |a2|  '); fprintf('  d=%-5.2f', dList); fprintf('\n');
for ir = 1:numel(relList)
    for ip = 1:numel(psiList)
        a2 = relList(ir)*exp(1j*psiList(ip)*pi/180);
        r = zeros(size(dList));
        for k = 1:numel(dList)
            g1 = gfun(x, 0); g2 = gfun(x, dList(k));
            J = [g1.*x/w^2, g1, 1j*g1, a2*g2.*(x - dList(k))/w^2, g2, 1j*g2];
            r(k) = infl(J);
        end
        if ir == 1, ratio2(ip, :) = r; end
        fprintf('  %6d   %4.1f ', psiList(ip), relList(ir)); fprintf('  %7.2f', r); fprintf('\n');
    end
end
dS = [0.05 0.1 0.2]; rS = zeros(2, 3);
for k = 1:3
    for ip = 1:2
        a2 = exp(1j*(ip - 1)*pi/2); g1 = gfun(x, 0); g2 = gfun(x, dS(k));
        rS(ip, k) = infl([g1.*x/w^2, g1, 1j*g1, a2*g2.*(x - dS(k))/w^2, g2, 1j*g2]);
    end
end
fprintf('  small-d slope of log(inflation) vs log(d): in phase %.2f, quadrature %.2f\n', ...
    log(rS(1,3)/rS(1,1))/log(4), log(rS(2,3)/rS(2,1))/log(4));

%% Part 3: carrier with known vs unknown phase (axial, z in wavelengths)
fwz = 0.75; wz = fwz/(2*sqrt(2*log(2))); z = (-3:1/32:3)'; k0 = 2*pi*2;   % two-way carrier: 2/lambda
g = exp(-z.^2/(2*wz^2)); s = g.*exp(1j*k0*z);
ds = (-z/wz^2 + 1j*k0).*s;                    % d s / d z0 (sign irrelevant for Fisher)
Iknown = 2*real(ds'*ds);                      % amplitude and phase known
Junk = [ds, s, 1j*s];                         % complex amplitude unknown
Ieff = 1/inv11(2*real(Junk'*Junk));
dg = -z/wz^2.*g; Ienv = 2*real(dg'*dg);       % envelope only, no carrier
fprintf('\nPart 3: axial Fisher (FWHM %.2f lambda, carrier 2/lambda)\n', fwz);
fprintf('  phase known: %.1f   phase unknown: %.2f   envelope only: %.2f   ratio unknown/envelope = %.4f\n', ...
    Iknown, Ieff, Ienv, Ieff/Ienv);
fprintf('  std gain from a known phase: %.1fx\n', sqrt(Iknown/Ieff));

%% Part 4: sound speed from receive delays of one point target
c = 1540; N = 128; pitch = 0.3e-3; sig = 1e-9; D = N*pitch;
xe = ((0:N-1)' - (N-1)/2)*pitch;
zList = [10 20 30 45 60 80]*1e-3;
sets = {'c only', 'c, z', 'c, t0', 'c, z, t0', 'c, z, xs, t0'};
sigc = zeros(numel(zList), numel(sets));
for iz = 1:numel(zList)
    zz = zList(iz); r = sqrt(zz^2 + xe.^2);
    dc = -(zz + r)/c^2; dz = (1 + zz./r)/c; dxs = xe./(r*c); dt0 = ones(N, 1);
    Js = {dc, [dc dz], [dc dt0], [dc dz dt0], [dc dz dxs dt0]};
    for k = 1:numel(Js)
        Ci = inv(Js{k}'*Js{k}/sig^2); sigc(iz, k) = sqrt(Ci(1, 1));
    end
end
fprintf('\nPart 4: sigma_c (m/s), %d elements, pitch %.1f mm, sigma_tau = %.1f ns\n', N, pitch*1e3, sig*1e9);
fprintf('  z(mm)  F#  '); fprintf('%14s', sets{:}); fprintf('   approx(c,z)\n');
for iz = 1:numel(zList)
    fprintf('  %4.0f  %4.2f', zList(iz)*1e3, zList(iz)/D); fprintf('%14.4f', sigc(iz, :));
    fprintf('   %9.4f\n', sqrt(180)*sig*zList(iz)*c^2/(sqrt(N)*D^2));
end

%% Figures
figure('Position', [80 80 950 360]);
subplot(1, 2, 1);
semilogy(dList, ratio2(1, :), 'b-o', dList, ratio2(2, :), 'r-s'); hold on; semilogy(dList, ones(size(dList)), 'k:');
xlabel('separation d (FWHM)'); ylabel('std inflation of x_1'); grid on;
title('Neighbour bubble, position and amplitude unknown'); legend('in phase', 'quadrature');
subplot(1, 2, 2);
semilogy(zList*1e3, sigc(:, 1), 'k-o', zList*1e3, sigc(:, 2), 'b-s', zList*1e3, sigc(:, 3), 'm-^', zList*1e3, sigc(:, 4), 'r-d');
xlabel('depth z (mm)'); ylabel('\sigma_c (m/s)'); grid on;
title('Sound speed from one target, 1 ns jitter'); legend(sets(1:4), 'Location', 'northwest');
