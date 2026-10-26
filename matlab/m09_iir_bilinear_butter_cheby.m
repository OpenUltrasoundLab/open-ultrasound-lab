%% m09_iir_bilinear_butter_cheby.m
% Module 09 - IIR design without a toolbox: analog prototype + prewarped bilinear
% transform, order formulas, biquads, and why direct form is fragile.
%  1) Orders for the RF anti-alias spec (fs = 40 MHz, pass 8 MHz with Ap = 0.17 dB,
%     stop 9 MHz with As = 60 dB): Butterworth, Chebyshev, elliptic (via ellipke)
%  2) Butterworth and Chebyshev I designs from their analog poles, bilinear with
%     prewarping, implemented as second-order sections; check the edges land exactly
%  3) Coefficient quantization: 8th-order Butterworth at fc = 0.02 fs in direct form
%     (16- and 24-bit words) vs biquads; 12th-order at 0.01 fs in double precision
% What to change: fp, fst, Ap, As; the word lengths; N and fc of the fragile filter.
% Expected output: orders 53 / 16 / 8; |H(8 MHz)| = -0.174 dB and |H(9 MHz)| <= -60 dB;
% direct-form max pole radius 1.299 (16 bit) and 1.166 (24 bit), biquads 0.9758;
% the N = 12 direct-form impulse response diverges, the biquad version decays.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.

clear; close all; clc;
rng(1);

%% Spec and order formulas
fs = 40; fp = 8; fst = 9;                  % MHz
Ap = 20*log10(1.01/0.99); As = 60;         % dB: the +-0.01 passband of the FIR spec, 60 dB stopband
ep = sqrt(10^(Ap/10) - 1); es = sqrt(10^(As/10) - 1);
Wp = 2*tan(pi*fp/fs); Ws = 2*tan(pi*fst/fs);          % prewarped edges (T = 1, 2/T = 2)
Nb = ceil(log(es/ep) / log(Ws/Wp));
Nc = ceil(acosh(es/ep) / acosh(Ws/Wp));
k  = Wp/Ws; k1 = ep/es;                                % ellipke takes m = k^2
[K, ~]  = ellipke(k^2);  [Kp, ~]  = ellipke(1 - k^2);
[K1, ~] = ellipke(k1^2); [K1p, ~] = ellipke(1 - k1^2);
Ne = ceil(K*K1p / (Kp*K1));
fprintf('Spec: pass %g MHz (Ap = %.3f dB), stop %g MHz (As = %g dB), fs = %g MHz\n', fp, Ap, fst, As, fs);
fprintf('  prewarped edges: Wp = %.4f, Ws = %.4f (selectivity %.4f)\n', Wp, Ws, Ws/Wp);
fprintf('  orders: Butterworth %d, Chebyshev I/II %d, elliptic %d\n', Nb, Nc, Ne);

%% Helpers (anonymous only)
bil   = @(s) (2 + s) ./ (2 - s);                          % bilinear map s -> z (2/T = 2)
evalz = @(zz, pp, g, w) g * prod(exp(1j*w(:).') - zz(:), 1) ./ prod(exp(1j*w(:).') - pp(:), 1);
wgrid = linspace(0, pi, 8001); f = wgrid/(2*pi)*fs;

%% Butterworth (passband edge exact) and Chebyshev I from analog poles
for fam = 1:2
    if fam == 1
        N = Nb; Wc = Wp / ep^(1/N); kk = 1:N;
        pa = Wc * exp(1j*pi*(2*kk + N - 1)/(2*N)); name = 'Butterworth';
    else
        N = Nc; mu = asinh(1/ep)/N; kk = 1:N; th = pi*(2*kk - 1)/(2*N);
        pa = Wp * (-sinh(mu)*sin(th) + 1j*cosh(mu)*cos(th)); name = 'Chebyshev I';
    end
    pd = bil(pa); zd = -ones(1, N);                        % N analog zeros at infinity -> z = -1
    dc = 1; if fam == 2 && mod(N, 2) == 0, dc = 1/sqrt(1 + ep^2); end
    g = real(dc * prod(1 - pd) / prod(1 - zd));            % set the DC gain
    H = abs(evalz(zd, pd, g, wgrid));
    % biquads: conjugate pole pairs with a double zero at z = -1
    pu = pd(imag(pd) > 1e-12); pr = pd(abs(imag(pd)) <= 1e-12);
    y = [1 zeros(1, 2999)]; gs = g;
    for i = 1:numel(pu)
        a2 = [1, -2*real(pu(i)), abs(pu(i))^2]; y = filter([1 2 1], a2, y);
    end
    for i = 1:numel(pr), y = filter([1 1], [1 -pr(i)], y); end
    y = gs * y;
    Hs = abs(exp(-1j*wgrid(:)*(0:numel(y)-1)) * y(:));     % response of the SOS impulse response
    fprintf('\n%s, N = %d (%d biquads):\n', name, N, ceil(N/2));
    fprintf('  |H| at %g MHz = %.3f dB (target -%.3f), at %g MHz = %.2f dB (target <= -%g)\n', fp, 20*log10(abs(evalz(zd, pd, g, 2*pi*fp/fs))), Ap, fst, 20*log10(abs(evalz(zd, pd, g, 2*pi*fst/fs))), As);
    fprintf('  max pole radius %.5f; SOS impulse response vs zpk: max |dH| = %.1e\n', max(abs(pd)), max(abs(Hs(:) - H(:))));
    if fam == 1, Hb = H; else, Hc = H; end
end

%% Coefficient quantization: direct form vs biquads
fprintf('\nCoefficient quantization (Butterworth low-pass, -3 dB at fc):\n');
N = 8; fc = 0.02; Wc = 2*tan(pi*fc); kk = 1:N;
pd = bil(Wc * exp(1j*pi*(2*kk + N - 1)/(2*N)));
a = real(poly(pd));
ib = ceil(log2(max(abs(a)))) + 1;                       % integer bits incl. sign
fprintf('  N = %d, fc = %.2f fs: true max pole radius %.4f, max |a_k| = %.1f (%d integer bits)\n', N, fc, max(abs(pd)), max(abs(a)), ib);
for W = [16 24]
    fr = W - ib; aq = round(a * 2^fr) / 2^fr;
    rq = max(abs(roots(aq)));
    pu = pd(imag(pd) > 0); r2 = 0;
    for i = 1:numel(pu)
        sec = [1, -2*real(pu(i)), abs(pu(i))^2]; secq = round(sec * 2^(W-2)) / 2^(W-2); r2 = max(r2, max(abs(roots(secq))));
    end
    fprintf('  %d-bit words: direct form max pole radius %.4f %s; biquads (%d fractional bits) %.4f\n', W, rq, char('(unstable)' .* (rq > 1) + ' ' .* (rq <= 1)), W-2, r2);
    if W == 16, pq16 = roots(aq); end
end
N = 12; fc = 0.01; Wc = 2*tan(pi*fc); kk = 1:N;
pd = bil(Wc * exp(1j*pi*(2*kk + N - 1)/(2*N))); a = real(poly(pd)); g = real(prod(1 - pd) / 2^N); b = g * real(poly(-ones(1, N)));
x = [1 zeros(1, 39999)]; yd = filter(b, a, x);
ys = x; pu = pd(imag(pd) > 0); for i = 1:numel(pu), ys = filter([1 2 1], [1, -2*real(pu(i)), abs(pu(i))^2], ys); end; ys = g * ys;
fprintf('  N = 12, fc = 0.01 fs, double precision: direct-form tail max |h| = %.3g (NaN/Inf = diverged); biquad tail max |h| = %.3g\n', max(abs(yd(end-1999:end))), max(abs(ys(end-1999:end))));

%% Figures
figure('Name', 'IIR by bilinear transform');
subplot(2,1,1); plot(f, 20*log10(Hb + 1e-15), 'k', f, 20*log10(Hc + 1e-15), 'b', 'LineWidth', 1.2); hold on;
plot([fst fs/2], [-As -As], 'r--', [fp fp], [-120 5], 'r:', [fst fst], [-120 5], 'r:');
ylim([-100 5]); xlim([0 20]); grid on; xlabel('frequency (MHz)'); ylabel('|H| (dB)'); title(sprintf('Butterworth N = %d (black), Chebyshev I N = %d (blue)', Nb, Nc));
subplot(2,1,2); N = 8; Wc = 2*tan(pi*0.02); kk = 1:N; pt = bil(Wc * exp(1j*pi*(2*kk + N - 1)/(2*N)));
tt = linspace(-0.6, 0.6, 200); plot(cos(tt), sin(tt), 'k-'); hold on;
plot(real(pt), imag(pt), 'bx', 'MarkerSize', 9, 'LineWidth', 1.5); plot(real(pq16), imag(pq16), 'rx', 'MarkerSize', 9, 'LineWidth', 1.5);
axis image; axis([0.4 1.45 -0.5 0.5]); grid on; xlabel('Re z'); ylabel('Im z'); title('N = 8 at 0.02 fs: true poles (blue), 16-bit direct form (red)');
