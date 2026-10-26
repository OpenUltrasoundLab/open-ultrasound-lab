%% m03_convolution_by_hand.m
% Purpose : compute one convolution four ways (input-side scatter loop,
%           output-side flip-and-slide loop, conv, convolution matrix),
%           show the wrap-around of a too-short FFT, check the convolution
%           theorem at one frequency, and time loop vs conv vs FFT.
% Change  : x, h, Nshort; Lr and the filter length in the timing cell.
% Expect  : all four methods print [2 5 7 1 -3]; the 4-point FFT gives
%           [-1 5 7 1]; X*H = Y = -8 - 4j at w = pi/2; sum(y) = 12 = 6*2.
clear; close all; rng(1);

x = [1 2 3];          % input
h = [2 1 -1];         % impulse response
L = numel(x); M = numel(h); Ny = L + M - 1;

%% 1. Input side (scatter): every x(k) launches a scaled, shifted copy of h
y1 = zeros(1, Ny);
for k = 1:L
    y1(k:k+M-1) = y1(k:k+M-1) + x(k) * h;
end

%% 2. Output side (gather): y(n) = sum_k x(k) h(n-k), h flipped and slid
y2 = zeros(1, Ny);
for n = 1:Ny
    for k = 1:L
        m = n - k + 1;                 % index into h (1-based)
        if m >= 1 && m <= M
            y2(n) = y2(n) + x(k) * h(m);
        end
    end
end

%% 3. Built-in conv, 4. convolution (Toeplitz) matrix
y3 = conv(x, h);
T = zeros(Ny, L);
for k = 1:L
    T(k:k+M-1, k) = h(:);              % column k = h shifted down by k-1
end
y4 = (T * x(:)).';
fprintf('scatter loop : %s\n', mat2str(y1));
fprintf('gather loop  : %s\n', mat2str(y2));
fprintf('conv         : %s\n', mat2str(y3));
fprintf('matrix T*x   : %s\n', mat2str(y4));
fprintf('length L+M-1 = %d, sum(y) = %g = sum(x)*sum(h) = %g\n', Ny, sum(y3), sum(x)*sum(h));
fprintf('alternating sum (w = pi): %g = %g\n', sum(y3 .* (-1).^(0:Ny-1)), ...
        sum(x .* (-1).^(0:L-1)) * sum(h .* (-1).^(0:M-1)));

%% Circular convolution: FFT products wrap when N < L+M-1
Nshort = 4;
yc4 = real(ifft(fft(x, Nshort) .* fft(h, Nshort)));
yc8 = real(ifft(fft(x, 8) .* fft(h, 8)));
fprintf('FFT, N = %d   : %s   (y(5) = -3 wrapped onto y(1))\n', Nshort, mat2str(round(yc4 * 1e10) / 1e10));
fprintf('FFT, N = 8   : %s\n', mat2str(round(yc8(1:Ny) * 1e10) / 1e10));

%% Convolution theorem at w = pi/2
dtft = @(s, w) sum(s .* exp(-1j * w * (0:numel(s)-1)));
w0 = pi/2;
X = dtft(x, w0); H = dtft(h, w0); Y = dtft(y3, w0);
fprintf('w = pi/2: X = %s, H = %s\n', num2str(X), num2str(H));
fprintf('          X*H = %s, Y from y = %s\n', num2str(X*H), num2str(Y));
fprintf('          |X||H| = %.3f = |Y| = %.3f; angle X + angle H = %.1f deg, angle Y = %.1f deg\n', ...
        abs(X)*abs(H), abs(Y), (angle(X)+angle(H))*180/pi, angle(Y)*180/pi);

%% Timing on an RF-sized problem (101-tap filter, 20000 samples)
Lr = 20000; Mf = 101;
r = randn(1, Lr); hb = randn(1, Mf);
tic; ya = zeros(1, Lr+Mf-1);
for m = 1:Mf                            % scatter over taps: Mf vector adds
    ya(m:m+Lr-1) = ya(m:m+Lr-1) + hb(m) * r;
end
t1 = toc;
tic; yb = conv(r, hb); t2 = toc;
Nf = 2^nextpow2(Lr+Mf-1);
tic; yc = real(ifft(fft(r, Nf) .* fft(hb, Nf))); yc = yc(1:Lr+Mf-1); t3 = toc;
fprintf('timing (s): tap loop %.4f, conv %.4f, FFT (N = %d) %.4f\n', t1, t2, Nf, t3);
fprintf('max |difference|: loop-conv %.1e, FFT-conv %.1e\n', max(abs(ya-yb)), max(abs(yc-yb)));

%% Figures
figure('Name', 'Convolution by hand');
cols = [0.80 0.55 0.10; 0.80 0.20 0.45; 0.10 0.55 0.40];
subplot(2,2,1); stem(0:L-1, x, 'filled'); hold on; stem(0:M-1, h, 'r'); hold off;
xlim([-0.5 4.5]); title('x (blue) and h (red)'); xlabel('n (samples)');
subplot(2,2,2); hold on;
for k = 1:L
    stem((k-1:k+M-2) + 0.12*(k-2), x(k)*h, 'Color', cols(k,:), 'filled');
end
hold off; xlim([-0.5 4.5]); title('copies x(k) h(n-k)'); xlabel('n (samples)');
subplot(2,2,3); stem(0:Ny-1, y3, 'k', 'filled'); hold on; stem(0:Nshort-1, yc4, 'r'); hold off;
xlim([-0.5 4.5]); title('linear (k), 4-pt circular (r)'); xlabel('n (samples)');
subplot(2,2,4); w = linspace(0, pi, 400);
Xw = arrayfun(@(v) dtft(x, v), w); Hw = arrayfun(@(v) dtft(h, v), w); Yw = arrayfun(@(v) dtft(y3, v), w);
plot(w/pi, abs(Xw), 'b', w/pi, abs(Hw), 'r', w/pi, abs(Yw), 'k', w/pi, abs(Xw.*Hw), 'g--');
xlabel('\omega/\pi'); ylabel('magnitude'); title('|X|,|H|,|Y|; |X||H| dashed');
