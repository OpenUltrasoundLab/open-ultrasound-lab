%% m09_wall_filters_ensemble.m
% Module 09 - Doppler wall (clutter) filters on a short slow-time ensemble, judged as
% matrices y = A x. Power response |H(f)|^2 = ||A e_f||^2 / N_out (Bjaerum et al. 2002).
%  FIR  : order K windowed-sinc high-pass (Hamming), N - K valid outputs
%  IIR  : 2nd-order Butterworth high-pass (bilinear), zero / step / projection init
%  Regr.: polynomial regression of degree p, A = I - Q Q'
% What to change: N, fc, K, p, the tissue Doppler ft and the clutter level.
% Expected output (N = 10, fc = 0.1 PRF): at DC the zero-state IIR gives only -12.6 dB;
% at 0.01 PRF step init -32 dB, projection -43 dB, regression p = 2 -63 dB;
% regression -3 dB cutoffs ~ (0.44 + 0.41 p) PRF/N.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.

clear; close all; clc;
rng(1);

N  = 10;          % ensemble length (pulses)
fc = 0.1;         % IIR / FIR cutoff, normalized to the PRF
K  = 4;           % FIR order
p  = 2;           % regression degree

%% Filters
% IIR: Butterworth low-pass prototype -> high-pass (s -> Wc/s) -> bilinear (2/T = 2)
no = 2; Wc = 2*tan(pi*fc); kk = 1:no;
pl = exp(1j*pi*(2*kk + no - 1)/(2*no));                   % prototype poles (Wc = 1)
ph = Wc ./ pl;                                             % high-pass poles; zeros at s = 0
pz = (2 + ph) ./ (2 - ph); a = real(poly(pz)); b = real(poly(ones(1, no)));
b = b * abs(sum(a .* (-1).^(0:no))) / abs(sum(b .* (-1).^(0:no)));   % unit gain at Nyquist
% steady-state initial state for a unit step (direct form II transposed, as filter() uses)
nz = numel(a) - 1; Mz = eye(nz); for i = 1:nz, Mz(i,1) = Mz(i,1) + a(i+1); if i < nz, Mz(i,i+1) = Mz(i,i+1) - 1; end; end
zi = Mz \ (b(2:end) - a(2:end)*b(1))';
% FIR: windowed-sinc high-pass, low-pass normalized to unit DC gain (exact zero at DC)
nf = (0:K)'; mf = nf - K/2; lp = 2*fc*ones(K+1, 1); q = mf ~= 0; lp(q) = sin(2*pi*fc*mf(q)) ./ (pi*mf(q));
lp = lp .* (0.54 - 0.46*cos(2*pi*nf/K)); hfir = -lp / sum(lp); hfir(K/2+1) = hfir(K/2+1) + 1;

%% Matrices
I = eye(N);
F = zeros(N); S = zeros(N); for k = 1:N, F(:,k) = filter(b, a, I(:,k)); S(:,k) = filter(b, a, I(:,k), zi*I(1,k)); end
G = zeros(N, nz); for i = 1:nz, e0 = zeros(nz, 1); e0(i) = 1; G(:,i) = filter(b, a, zeros(N,1), e0); end
P = (I - G*pinv(G)) * F;                                    % projection initialization
t = linspace(-1, 1, N)'; V = zeros(N, p+1); for d = 0:p, V(:,d+1) = t.^d; end; [Q, ~] = qr(V, 0);
R = I - Q*Q';                                               % polynomial regression
Af = zeros(N-K, N); for i = 1:N-K, Af(i, i:i+K) = fliplr(hfir'); end   % FIR, valid outputs only

presp = @(A, f) sum(abs(A * exp(2j*pi*(0:N-1)'*f(:)')).^2, 1) / size(A, 1);
fr = [0 0.01 0.02 0.05 0.25];
names = {'FIR K=4', 'IIR zero init', 'IIR step init', 'IIR projection', 'regression p=2'};
mats = {Af, F, S, P, R};
fprintf('Power response (dB), N = %d, fc = %.2f PRF\n', N, fc);
fprintf('%-16s %6s %7s %7s %7s %7s   Nout\n', 'filter', 'f=0', '0.01', '0.02', '0.05', '0.25');
for i = 1:numel(mats)
    r = 10*log10(presp(mats{i}, fr) + 1e-30); r(r < -99) = -99;
    fprintf('%-16s %6.1f %7.1f %7.1f %7.1f %7.1f   %d\n', names{i}, r, size(mats{i}, 1));
end

%% Regression cutoff versus degree and ensemble length
ff = linspace(0, 0.5, 20001);
fprintf('\nRegression -3 dB cutoff x N / PRF  (approx. 0.44 + 0.41 p):\n     N   p=0    p=1    p=2    p=3\n');
for NN = [8 10 12 16 32]
    tt = linspace(-1, 1, NN)'; row = zeros(1, 4);
    for pp = 0:3
        Vv = zeros(NN, pp+1); for d = 0:pp, Vv(:,d+1) = tt.^d; end; [Qq, ~] = qr(Vv, 0);
        resp = 1 - sum(abs(Qq' * exp(2j*pi*(0:NN-1)'*ff)).^2, 1) / NN;
        row(pp+1) = ff(find(resp >= 0.5, 1)) * NN;
    end
    fprintf('  %4d  %5.3f  %5.3f  %5.3f  %5.3f\n', NN, row);
end
fprintf('  N = 10, p = 2 at PRF 4 kHz, f0 = 5 MHz: cutoff %.0f Hz = %.1f cm/s\n', 1.267*4000/10, 100*1540*1.267*4000/10/(2*5e6));

%% Clutter + blood ensemble
ft = 0.01; n = (0:N-1)';
c = 100 * exp(2j*pi*ft*n);                                  % clutter, 40 dB above blood
bl = exp(2j*pi*0.25*n + 0.7j);                              % blood at 0.25 PRF
fprintf('\nClutter-to-blood ratio after filtering (input +40.0 dB), tissue at %.2f PRF:\n', ft);
for i = 1:numel(mats)
    A = mats{i}; fprintf('  %-16s %6.1f dB\n', names{i}, 10*log10(sum(abs(A*c).^2) / sum(abs(A*bl).^2)));
end

%% Figures
figure('Name', 'Wall filters');
subplot(2,1,1); fg = linspace(0, 0.5, 400); cols = [0.7 0.45 0; 0.5 0.5 0.5; 0.85 0.65 0.1; 0.8 0.2 0.4; 0.1 0.35 0.8];
hold on; for i = 1:numel(mats), plot(fg, 10*log10(presp(mats{i}, fg) + 1e-12), 'Color', cols(i,:), 'LineWidth', 1.3); end
ylim([-80 5]); grid on; xlabel('f / PRF'); ylabel('power response (dB)'); title(sprintf('N = %d, f_c = %.2f PRF', N, fc)); legend(names, 'Location', 'southeast');
subplot(2,1,2); hold on; x = c + bl;
for i = 1:numel(mats)
    y = mats{i} * x; off = N - size(mats{i}, 1); plot((0:numel(y)-1) + off, 20*log10(abs(y)), 'o-', 'Color', cols(i,:), 'LineWidth', 1.2);
end
plot([0 N-1], [0 0], 'k:'); ylim([-40 45]); grid on; xlabel('slow-time sample n'); ylabel('|y[n]| (dB)'); title('Clutter (40 dB, 0.01 PRF) + blood (0 dB, 0.25 PRF)');
