%% m02_fft_radix2.m  -- Module 02: radix-2 FFT, recursive and iterative, checked against fft
% Purpose : implement X[k] = E[k] + W^k O[k], X[k+N/2] = E[k] - W^k O[k] two ways:
%           (a) recursive, written with anonymous functions only (no local functions), and
%           (b) iterative in place: bit-reversal permutation, then log2(N) butterfly stages.
%           Check both against fft, show the bit-reversed order and time them against the
%           direct O(N^2) DFT.
% Change  : Nlist (powers of two), the test signal.
% Expect  : bit-reversed order for N = 8: 0 4 2 6 1 5 3 7; N = 4 by hand E = [4 -2], O = [6 -2],
%           X = [10, -2+2j, -2, -2-2j]; errors vs fft ~1e-13 or smaller; direct DFT time grows
%           ~4x per doubling of N, FFTs ~2x; multiplies (N/2)log2(N) = 5120 vs N^2 = 1048576
%           at N = 1024 (ratio 204.8).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

%% 1. Recursive radix-2 FFT with anonymous functions
% A cell array of two handles plays the role of if/else: only the selected branch is called.
tw     = @(N) exp(-2j*pi*(0:N/2-1).'/N);                      % twiddle factors W_N^k, k = 0..N/2-1
comb   = @(E, O, N) [E + tw(N).*O; E - tw(N).*O];             % one level of butterflies
branch = {@(f, x) x, @(f, x) comb(f(f, x(1:2:end)), f(f, x(2:2:end)), numel(x))};
rfft   = @(f, x) branch{(numel(x) > 1) + 1}(f, x);            % 1-point DFT is the sample itself
fft_rec = @(x) rfft(rfft, x(:));

x4 = [1; 2; 3; 4];
E = fft([1; 3]); O = fft([2; 4]);
fprintf('--- N = 4 by hand ---\n');
fprintf('E = DFT[1 3] = [%g %g], O = DFT[2 4] = [%g %g], W4^1 = -j\n', E, O);
X4 = fft_rec(x4);
fprintf('recursive FFT: X = [%s]\n', num2str(X4.', '%6.2f'));

%% 2. Iterative in-place radix-2 FFT
bitrev = @(N) bin2dec(fliplr(dec2bin(0:N-1, round(log2(N))))) + 1;   % 1-based bit-reversed index
fprintf('\n--- bit reversal ---\nN = 8 input order: %s\n', num2str((bitrev(8) - 1).'));
fprintf('N = 16 input order: %s\n', num2str((bitrev(16) - 1).'));

%% 3. Accuracy against fft
Nlist = 2.^(6:12);
nN = numel(Nlist);
t_dft = nan(1, nN); t_rec = zeros(1, nN); t_it = zeros(1, nN); t_fft = zeros(1, nN);
err_rec = zeros(1, nN); err_it = zeros(1, nN);
fprintf('\n     N   err(recursive)  err(iterative)  t_dft(ms)  t_rec(ms)  t_iter(ms)  t_fft(ms)\n');
for i = 1:nN
    N = Nlist(i); L = round(log2(N));
    x = randn(N, 1) + 1j*randn(N, 1);
    Xref = fft(x);
    % recursive
    tic; Xr = fft_rec(x); t_rec(i) = toc;
    % iterative: permute once, then stages; group the butterflies of each stage as columns
    tic;
    X = x(bitrev(N));
    for s = 1:L
        m = 2^s; h = m/2;
        w = exp(-2j*pi*(0:h-1).'/m);
        Xs = reshape(X, m, N/m);                     % each column is one group of length m
        a = Xs(1:h, :); b = w .* Xs(h+1:m, :);       % implicit expansion of w over the groups
        Xs = [a + b; a - b];
        X = Xs(:);
    end
    t_it(i) = toc;
    % direct DFT by matrix (only up to N = 2048 to keep memory small)
    if N <= 2048
        tic; [kk, nn] = ndgrid(0:N-1, 0:N-1); Xd = exp(-2j*pi*kk.*nn/N)*x; t_dft(i) = toc; %#ok<NASGU>
    end
    tic; for r = 1:20, Xf = fft(x); end; t_fft(i) = toc/20; %#ok<NASGU>
    err_rec(i) = max(abs(Xr - Xref))/max(abs(Xref));
    err_it(i)  = max(abs(X  - Xref))/max(abs(Xref));
    fprintf('%6d   %12.1e   %12.1e   %9.2f  %9.2f  %10.2f  %9.3f\n', N, err_rec(i), err_it(i), ...
            1e3*t_dft(i), 1e3*t_rec(i), 1e3*t_it(i), 1e3*t_fft(i));
end

%% 4. Operation counts
fprintf('\n     N        N^2   (N/2)log2N   ratio 2N/log2N\n');
for N = [8 256 1024 4096 65536]
    fprintf('%6d %10d %12d %12.1f\n', N, N^2, N/2*log2(N), N^2/(N/2*log2(N)));
end

%% 5. Figure: timing and counts
figure('Position', [100 100 900 380]);
subplot(1, 2, 1);
loglog(Nlist, 1e3*t_dft, 'o-', Nlist, 1e3*t_rec, 's-', Nlist, 1e3*t_it, 'd-', Nlist, 1e3*t_fft, '^-');
xlabel('N'); ylabel('time (ms)'); grid on; legend('direct DFT (matrix)', 'recursive radix-2', 'iterative radix-2', 'built-in fft', 'Location', 'northwest');
title('Run time (interpreter overhead dominates the hand-written FFTs)');
subplot(1, 2, 2);
Nc = 2.^(2:16);
loglog(Nc, Nc.^2, 'o-', Nc, Nc/2.*log2(Nc), 's-'); grid on;
xlabel('N'); ylabel('complex multiplies'); legend('N^2', '(N/2) log_2 N', 'Location', 'northwest'); title('Operation count');
