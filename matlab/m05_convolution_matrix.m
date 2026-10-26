%% m05_convolution_matrix.m  Module 05: convolution is a Toeplitz matrix, the DFT diagonalizes circulants
% Purpose : build the Toeplitz convolution matrix T of a Gaussian-modulated pulse h,
%           verify T*x == conv(h,x); embed it in a circulant C, show that the DFT
%           matrix diagonalizes C with fft(first column) on the diagonal, and show
%           the circular wrap-around error when the circulant is too small.
% Change  : f0, frac_bw, fs (pulse); N (signal length); Nc (circulant size).
% Expect  : square 6x6 window: T6*x6 = first 6 samples of conv, last 3 lost;
%           max|T*x - conv(h,x)| = 0; max off-diagonal |F C F'/Nc| ~ 1e-15;
%           with Nc = N the circular result differs from linear conv in the first L-1 samples;
%           direct vs FFT convolution of a 2048-sample line cross over near L ~ 100 taps.
% Runtime : ~1 s.  Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);

%% 1. A square 6x6 convolution window (h has 4 taps): the full result needs 6+4-1 = 9 rows
h4 = [-0.23; -0.62; 0.62; 0.23];         % sampled Gaussian-modulated sine, zero mean
x6 = [1; 0; 0; 0.5; 0; -1];               % sparse reflectivity, 6 samples
T6 = zeros(6,6);
for j = 1:6                               % column j = h shifted down by j-1 (truncated at 6 rows)
    rows = j:min(6, j+3);
    T6(rows, j) = h4(1:numel(rows));
end
disp('T6 (lower-triangular Toeplitz, T(i,j) = h(i-j+1)) ='); disp(T6);
yfull = conv(h4, x6);                     % full linear convolution: 6 + 4 - 1 = 9 samples
fprintf('T6*x6         = [%s]\n', num2str((T6*x6).', '%7.3f'));
fprintf('conv(h,x)(1:6)= [%s]\n', num2str(yfull(1:6).', '%7.3f'));
fprintf('tail conv(7:9)= [%s]  (lost by the square 6x6 window)\n', num2str(yfull(7:9).', '%7.3f'));

%% 2. A realistic pulse: f0 = 5 MHz, 60 percent -6 dB bandwidth, fs = 40 MHz
f0 = 5e6; frac_bw = 0.6; fs = 40e6;       % Hz, fraction, Hz
sig_f = frac_bw*f0/2 / sqrt(2*log(2));    % std of the Gaussian amplitude spectrum (Hz)
sig_t = 1/(2*pi*sig_f);                   % std of the Gaussian envelope (s)
t = (-ceil(3*sig_t*fs):ceil(3*sig_t*fs)).' / fs;   % +-3 sigma_t
h = exp(-t.^2/(2*sig_t^2)) .* cos(2*pi*f0*t);     % RF pulse
L = numel(h);
fprintf('pulse: sigma_t = %.4f us, L = %d taps at %.0f MHz\n', sig_t*1e6, L, fs/1e6);

%% 3. Full Toeplitz convolution matrix, (N+L-1) x N
N = 64;                                   % signal length (samples)
x = zeros(N,1); x([10 25 26 47]) = [1 -0.6 0.8 0.5];   % four scatterers
T = zeros(N+L-1, N);
for j = 1:N
    T(j:j+L-1, j) = h;                    % every column is the pulse, shifted by one sample
end
fprintf('T is %d x %d, %d nonzeros (%.1f %% dense), Toeplitz: T(i,j) depends only on i-j\n', ...
        size(T,1), size(T,2), nnz(T), 100*nnz(T)/numel(T));
fprintf('max|T*x - conv(h,x)| = %.1e\n', max(abs(T*x - conv(h,x))));

%% 4. Circulant embedding: zero-pad to Nc >= N+L-1 and linear = circular convolution
Nc = N + L - 1;
cfirst = [h; zeros(Nc-L,1)];              % first column of the circulant
C = zeros(Nc);
for j = 1:Nc
    C(:,j) = circshift(cfirst, j-1);      % each column = previous column rotated down by one
end
xp = [x; zeros(Nc-N,1)];
yc = real(ifft(fft(cfirst) .* fft(xp)));  % FFT convolution
fprintf('Nc = %d: max|C*xp - conv| = %.1e, max|ifft(fft.*fft) - conv| = %.1e\n', ...
        Nc, max(abs(C*xp - conv(h,x))), max(abs(yc - conv(h,x))));

%% 5. The DFT matrix diagonalizes the circulant
F = fft(eye(Nc));                         % DFT matrix: F(k+1,n+1) = exp(-2i*pi*k*n/Nc)
fprintf('F''*F = Nc*I: max|F''*F/Nc - I| = %.1e\n', max(max(abs(F'*F/Nc - eye(Nc)))));
Dg = F * C * F' / Nc;                     % should be diag(fft(cfirst))
offd = Dg - diag(diag(Dg));
fprintf('max off-diagonal |F C F''/Nc| = %.1e, max|diag - fft(c)| = %.1e\n', ...
        max(abs(offd(:))), max(abs(diag(Dg) - fft(cfirst))));
fprintf('eigenvalues of C = fft(c): |lambda| max %.3f at %.2f MHz\n', ...
        max(abs(fft(cfirst))), (find(abs(fft(cfirst)) == max(abs(fft(cfirst))), 1) - 1)*fs/Nc/1e6);

%% 6. Too small a circulant: circular wrap-around (time aliasing)
Ns = N;                                   % no zero padding
cs = [h; zeros(Ns-L,1)];
ycirc = real(ifft(fft(cs) .* fft(x)));
ylin = conv(h, x);
err = abs(ycirc - ylin(1:Ns));
fprintf('Nc = N = %d: max wrap error = %.3f in samples 1..%d, %.1e afterwards\n', ...
        Ns, max(err(1:L-1)), L-1, max(err(L:end)));

%% 7. Cost of the two routes for one RF line
Nrf = 2048; Lf = [16 64 256];
for k = 1:3
    Nfft = 2^nextpow2(Nrf + Lf(k) - 1);
    direct = Nrf * Lf(k);                                 % real multiply-adds
    viafft = 2 * 2*Nfft*log2(Nfft) + 4*Nfft;              % 2 real-equivalent FFTs (~2 N log2 N real mults each) + products
    fprintf('L = %3d: direct %7d MACs, FFT route ~%7d mults (Nfft = %d)\n', Lf(k), direct, round(viafft), Nfft);
end

%% 8. Figures
figure;
subplot(2,2,1); imagesc(T); colorbar; title('Toeplitz T (N+L-1 x N)'); xlabel('input sample j'); ylabel('output sample i');
subplot(2,2,2); imagesc(C); axis image; colorbar; title('circulant C (wraps at the corner)'); xlabel('column'); ylabel('row');
subplot(2,2,3); imagesc(20*log10(abs(Dg)/max(abs(Dg(:))) + 1e-12)); axis image; colorbar;
caxis([-300 0]); title('|F C F^H / N_c| (dB): diagonal'); xlabel('column k'); ylabel('row k');
subplot(2,2,4); plot(0:numel(ylin)-1, ylin, 'k', 0:Ns-1, ycirc, 'r--');
xlabel('sample'); ylabel('amplitude'); title('linear (black) vs circular N_c = N (red)');
colormap(gray);
