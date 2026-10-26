%% m23_broadband_beams.m  --  pulsed beams: grating-lobe level ~ M/N, far sidelobes = edge waves
% Purpose : broadband (pulsed) beam patterns of a uniform linear array, synthesized in the
%           frequency domain: for each direction u = sin(theta) the received pulse is
%           s(t) = IFFT{ P(f) AF(f,u)^m }, m = 1 (one-way) or 2 (two-way, Tx = Rx), with a
%           Gaussian pulse spectrum P(f). From s(t) we take the peak |s| (peak-pressure
%           pattern) and the energy sum |s|^2 (time-integrated pattern), and compare with:
%           (1) grating-lobe peak ~ 1.064 M/N (one- and two-way) and energy ~ 1.505 M/N
%               (one-way) or (2/3) 1.505 M/N = 1.003 M/N (two-way, triangular coarray), where
%               M = f0 x (-6 dB pulse duration) = 0.883/FBW cycles and N = element count;
%           (2) far sidelobes, peak ~ lambda/(2 pi D |u|): two separate edge waves, 6 dB
%               below the CW sidelobe envelope lambda/(pi D |u|).
% Change  : N, FBW, pitch (d/lambda), the Gaussian pulse model.
% Expect  : grating-lobe peaks within ~0.6 dB of 20 log10(1.064 M/N) for one-way and
%           two-way (N = 128, FBW 60 %: -38 dB, vs 0 dB in CW); energy within 0.2 dB of
%           its law. Far sidelobes within ~1 dB of the edge-wave law.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

c = 1.540; f0 = 5; lam = c/f0;          % mm/us, MHz, mm
fs = 40; nf = 4096; df = fs/nf;         % baseband frequency grid (MHz), time window nf/fs = 102 us
fb = (-nf/2:nf/2-1).' * df;             % baseband frequency f - f0 (MHz), column
% Uniform array factor in closed form (Dirichlet kernel), elements centered on 0
% psi = pi f d u / c ; AF = sin(N psi)/sin(psi)
sm     = @(psi) abs(sin(psi)) < 1e-9;                  % psi at a multiple of pi
dirich = @(psi, N) (~sm(psi)).*sin(N*psi)./(sin(psi) + sm(psi)) + sm(psi).*N.*cos(N*psi)./cos(psi);
fprintf('Gaussian pulse: M = 0.883/FBW cycles (-6 dB duration x f0)\n\n');

%% (1) Grating-lobe level: peak and energy, one-way and two-way
d = 1.5*lam;  ug = lam/d;               % first grating lobe at u = 2/3 (41.8 deg) at f0
fprintf('(1) grating lobe at u = %.3f (d = 1.5 lambda), level re main lobe (dB)\n', ug);
fprintf('%5s %5s %6s %5s | %8s %8s | %8s %8s\n', 'N', 'FBW', 'M', 'way', 'peak', '1.064M/N', 'energy', 'theory');
for N = [64 128]
    for FBW = [0.3 0.6 0.9]
        M = 0.883/FBW;
        for m = 1:2
            sf = FBW*f0/(2*sqrt(2*log(2)));            % std of the Gaussian amplitude spectrum
            P  = exp(-fb.^2/(2*sf^2));
            f  = f0 + fb;
            S  = zeros(nf, 2);
            uu = [0 ug];
            for j = 1:2
                AF = dirich(pi*f*d*uu(j)/c, N);
                S(:,j) = ifft(ifftshift(P .* AF.^m));
            end
            pk = max(abs(S));  en = sum(abs(S).^2);
            fprintf('%5d %5.1f %6.2f %5d | %8.1f %8.1f | %8.1f %8.1f\n', N, FBW, M, m, ...
                20*log10(pk(2)/pk(1)), 20*log10(1.064*M/N), 10*log10(en(2)/en(1)), 10*log10(1.505*M/N*(1 - (m == 2)/3)));
        end
    end
end

%% (2) Far sidelobes: edge waves (one-way, d = lambda/2, N = 128, FBW 60 %)
N = 128; d = lam/2; D = N*d; FBW = 0.6; M = 0.883/FBW;
sf = FBW*f0/(2*sqrt(2*log(2))); P = exp(-fb.^2/(2*sf^2)); f = f0 + fb;
uu = [0 0.05 0.1 0.2 0.4];
S  = zeros(nf, numel(uu));
for j = 1:numel(uu), S(:,j) = ifft(ifftshift(P .* dirich(pi*f*d*uu(j)/c, N))); end
pk = max(abs(S)); pk = pk/pk(1);
fprintf('\n(2) far sidelobes, one-way, N = %d, d = lambda/2 (D = %.1f mm), FBW %.0f %%\n', N, D, 100*FBW);
fprintf('    edge waves separate in time when |u| > M lambda/D = %.3f\n', M*lam/D);
fprintf('%7s %10s %14s %14s\n', 'u', 'peak (dB)', 'edge 1/(2x)', 'CW env 1/x');
for j = 2:numel(uu)
    x = pi*D*uu(j)/lam;
    fprintf('%7.2f %10.1f %14.1f %14.1f\n', uu(j), 20*log10(pk(j)), 20*log10(1/(2*x)), 20*log10(1/x));
end

%% (3) Patterns for the figure: N = 64, d = 1.5 lambda, FBW 60 %
N = 64; d = 1.5*lam; D = N*d; FBW = 0.6; M = 0.883/FBW;
sf = FBW*f0/(2*sqrt(2*log(2))); P = exp(-fb.^2/(2*sf^2)); f = f0 + fb;
u  = linspace(0, 1, 501);
pat = zeros(numel(u), 3, 2);               % [CW, peak, energy] x [one-way, two-way]
for m = 1:2
    for blk = 0:4
        jj = blk*101 + (1:101); jj = jj(jj <= numel(u));
        AF = dirich(pi*f*d*u(jj)/c, N);                     % nf x numel(jj)
        S  = ifft(ifftshift(P .* AF.^m, 1), [], 1);
        pat(jj, 2, m) = max(abs(S)).';
        pat(jj, 3, m) = sqrt(sum(abs(S).^2)).';
        pat(jj, 1, m) = abs(dirich(pi*f0*d*u(jj)/c, N)).'.^m;
    end
    for q = 1:3, pat(:,q,m) = pat(:,q,m)/pat(1,q,m); end
end
lab = {'one-way', 'two-way'};
figure('Position', [60 60 1150 440]);
for m = 1:2
    subplot(1,2,m);
    plot(u, 20*log10(pat(:,1,m) + 1e-9), 'Color', [0.6 0.6 0.6]); hold on;
    plot(u, 20*log10(pat(:,2,m) + 1e-9), 'LineWidth', 1.4);
    plot(u, 20*log10(pat(:,3,m) + 1e-9), '--', 'LineWidth', 1.2);
    plot(lam/d + [-0.06 0.06], [1 1]*20*log10(1.064*M/N), 'k-', 'LineWidth', 2);
    if m == 1, ue = u(u > 0.03); plot(ue, 20*log10(lam./(2*pi*D*ue)), 'k:', 'LineWidth', 1.2); end
    grid on; ylim([-80 3]); xlabel('u = sin \theta'); ylabel('dB re main lobe');
    title(sprintf('%s, N = %d, d = 1.5\\lambda, FBW %.0f %%', lab{m}, N, 100*FBW));
    if m == 1, legend('CW (f_0)', 'pulsed: peak', 'pulsed: energy (rms)', '1.064 M/N', 'edge waves \lambda/(2\piDu)', 'Location', 'northeast');
    else, legend('CW (f_0)', 'pulsed: peak', 'pulsed: energy (rms)', '1.064 M/N', 'Location', 'northeast'); end
end
