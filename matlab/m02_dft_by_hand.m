%% m02_dft_by_hand.m  -- Module 02: the DFT by hand for N = 4 and N = 8
% Purpose : compute X[k] = sum_n x[n] exp(-j*2*pi*k*n/N) three ways (explicit loops with every
%           multiply printed, the DFT matrix F, and fft), then read each bin as cosine content
%           (real part), minus sine content (imaginary part), amplitude and phase.
% Change  : x4 (the N = 4 example) and x8 (the N = 8 example).
% Expect  : x4 = [1 2 3 4] -> X = [10, -2+2j, -2, -2-2j]; x[n] = 2.5 - cos(pi n/2) - sin(pi n/2)
%           - 0.5(-1)^n; Parseval 30 = 120/4. Pulse [1 1 1 1 0 0 0 0] -> X[1] = 1 - 2.414j,
%           |X[1]| = 2.613, angle -67.5 deg, X[2] = X[4] = X[6] = 0. Loops, matrix and fft agree
%           to about 1e-15.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);
cfmt = @(z) sprintf('%7.3f %s %6.3fj', real(z), char(43 + 2*(imag(z) < 0)), abs(imag(z)));  % '+' or '-'

%% 1. N = 4: every multiply
x4 = [1 2 3 4];
N = numel(x4);
X_loop = zeros(1, N);
fprintf('--- N = 4, x = [%s] ---\n', num2str(x4));
fprintf('W = exp(-j*2*pi/4) = -j, so every probe value is 1, -j, -1 or j\n');
for k = 0:N-1
    terms = '';
    for n = 0:N-1
        W = exp(-1j*2*pi*mod(k*n, N)/N);                   % probe value W^(kn)
        W = round(W*1e12)/1e12;                            % clean -0 and 1e-16 for printing
        p = x4(n+1)*W;                                     % one multiply
        X_loop(k+1) = X_loop(k+1) + p;
        terms = [terms sprintf('  %d*(%s)', x4(n+1), strtrim(cfmt(W)))]; %#ok<AGROW>
    end
    fprintf('k = %d:%s\n        X[%d] = %s\n', k, terms, k, cfmt(X_loop(k+1)));
end

%% 2. The DFT matrix and fft
[kk, nn] = ndgrid(0:N-1, 0:N-1);
F4 = exp(-1j*2*pi*kk.*nn/N);                               % F(k+1, n+1) = W^(kn)
X_mat = (F4*x4(:)).';
X_fft = fft(x4);
fprintf('max |loop - fft| = %.1e, max |matrix - fft| = %.1e\n', max(abs(X_loop - X_fft)), max(abs(X_mat - X_fft)));
fprintf('F4''*F4 = %g * I (max off-diagonal %.1e)\n', real(trace(F4'*F4))/N, max(max(abs(F4'*F4 - N*eye(N)))));

%% 3. What each bin means
fprintf('\n--- meaning of each bin (N = 4) ---\n');
fprintf('X[0]/N = %.3f   (mean, DC)\n', real(X_fft(1))/N);
fprintf('bin 1: cos amplitude 2*Re/N = %.3f, sin amplitude -2*Im/N = %.3f, |X| = %.3f, angle = %.1f deg\n', ...
        2*real(X_fft(2))/N, -2*imag(X_fft(2))/N, abs(X_fft(2)), angle(X_fft(2))*180/pi);
fprintf('X[2]/N = %.3f   (Nyquist, multiplies (-1)^n; always real for real x)\n', real(X_fft(3))/N);
fprintf('X[3] = %s = conj(X[1]) (negative frequency)\n', cfmt(X_fft(4)));
n = 0:N-1;
x_rebuilt = real(X_fft(1))/N + (2/N)*(real(X_fft(2))*cos(2*pi*n/N) - imag(X_fft(2))*sin(2*pi*n/N)) + real(X_fft(3))/N*(-1).^n;
fprintf('rebuilt from the pieces: [%s]\n', num2str(x_rebuilt, '%6.3f'));
fprintf('Parseval: sum|x|^2 = %g, sum|X|^2/N = %g\n', sum(abs(x4).^2), sum(abs(X_fft).^2)/N);

%% 4. N = 8: the 4-sample pulse, table for bin 1
x8 = [1 1 1 1 0 0 0 0];
N = numel(x8); k = 1;
fprintf('\n--- N = 8, x = [%s], bin k = %d ---\n', num2str(x8), k);
fprintf('  n  x[n]  theta(deg)   cos     -sin    x*cos   -x*sin   sum Re   sum Im\n');
sre = 0; sim = 0;
for n = 0:N-1
    th = 2*pi*mod(k*n, N)/N;
    sre = sre + x8(n+1)*cos(th); sim = sim - x8(n+1)*sin(th);
    fprintf('  %d  %4.1f  %6.0f   %7.3f %7.3f %7.3f %7.3f  %7.3f  %7.3f\n', n, x8(n+1), th*180/pi, ...
            cos(th), -sin(th), x8(n+1)*cos(th), -x8(n+1)*sin(th), sre, sim);
end
X8 = fft(x8);
fprintf('X[1] = %s, |X[1]| = %.3f, angle = %.1f deg (pulse centred at n = 1.5: -2*pi*1.5/8 = %.1f deg)\n', ...
        cfmt(X8(2)), abs(X8(2)), angle(X8(2))*180/pi, -2*pi*1.5/8*180/pi);
fprintf('all 8 bins:\n');
for kk8 = 0:N-1
    fprintf('  X[%d] = %s   |X| = %.4f\n', kk8, cfmt(X8(kk8+1)), abs(X8(kk8+1)));
end
[kk, nn] = ndgrid(0:N-1, 0:N-1); F8 = exp(-1j*2*pi*kk.*nn/N);
fprintf('max |F8*x - fft(x)| = %.1e\n', max(abs((F8*x8(:)).' - X8)));

%% 5. Figures
figure('Position', [100 100 900 600]);
subplot(2, 2, 1); stem(0:7, x8, 'filled'); xlabel('sample n'); ylabel('x[n]'); title('4-sample pulse, N = 8'); xlim([-0.5 7.5]);
subplot(2, 2, 2); stem((0:7) - 0.1, real(X8), 'filled'); hold on; stem((0:7) + 0.1, imag(X8), 'filled'); hold off;
xlabel('bin k'); ylabel('X[k]'); legend('Re X (cosine content)', 'Im X (-sine content)'); title('Re even, Im odd about k = 4'); xlim([-0.5 7.5]);
subplot(2, 2, 3); imagesc(0:7, 0:7, real(F8)); axis image; colorbar; xlabel('n'); ylabel('k'); title('Re F_8 = cos(2\pi kn/8)');
subplot(2, 2, 4); imagesc(0:7, 0:7, imag(F8)); axis image; colorbar; xlabel('n'); ylabel('k'); title('Im F_8 = -sin(2\pi kn/8)');
