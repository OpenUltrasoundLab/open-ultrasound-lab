%% m05_das_as_matrix.m  Module 05: delay-and-sum as a sparse matrix, and its adjoint
% Purpose : build the DAS operator B (pixels x channel samples) for a 32-element
%           linear array and a 0-degree plane-wave transmit, with linear-interpolation
%           weights (2 nonzeros per element per pixel); apply it to simulated
%           point-target channel data; show that B' (the adjoint) maps an image back to
%           channel data along hyperbolas; check <B d, g> = <d, B' g>; spy(B) and density.
% Change  : Ne, pitch, the pixel grid (dxp), apod ('rect' or 'hann'), the targets.
% Expect  : B is 6161 x 43520 with 64 nonzeros per row (0.147 percent dense);
%           dot-test mismatch ~1e-15; lateral -6 dB width at z = 18 mm close to
%           1.21*lambda*z/D = 0.70 mm (rect); pulse-filtered B' column vs true echo
%           correlation ~0.99 (B is H^H up to the pulse matched filter).
% Runtime : ~3 s.  Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);

%% 1. System and grid
c = 1540;  f0 = 5e6;  fs = 40e6;  frac_bw = 0.6;    % m/s, Hz, Hz, -6 dB fractional bandwidth
lambda = c/f0;                                      % 0.308 mm
Ne = 32;   pitch = 0.3e-3;                          % elements, m
xe = ((0:Ne-1) - (Ne-1)/2) * pitch;                 % element x positions (m), aperture D = 9.6 mm
D = Ne*pitch;
apod = 'rect';
if strcmp(apod, 'hann'), we = sin(pi*((0:Ne-1) + 0.5)/Ne).^2; else, we = ones(1, Ne); end
dxp = 0.1e-3;                                       % pixel size (m), about lambda/3
xg = -3e-3:dxp:3e-3;  zg = 15e-3:dxp:25e-3;         % 61 x 101 pixels
[Xp, Zp] = meshgrid(xg, zg);  Np = numel(Xp);
Nt = 1360;                                          % samples per channel (34 us at 40 MHz)
t = (0:Nt-1).'/fs;
fprintf('lambda = %.3f mm, D = %.1f mm, pixels = %d x %d = %d, channel samples = %d x %d = %d\n', ...
        lambda*1e3, D*1e3, numel(zg), numel(xg), Np, Nt, Ne, Nt*Ne);

%% 2. Build the sparse DAS matrix B: row p = pixel p, column (e-1)*Nt + n = sample n of element e
Ip = zeros(2*Ne*Np, 1); Jc = Ip; Vw = Ip; k = 0;
for e = 1:Ne
    tau = (Zp(:) + sqrt((Xp(:) - xe(e)).^2 + Zp(:).^2)) / c;   % two-way delay, 0-degree plane wave (s)
    u = tau*fs;                                                 % fractional sample index (0-based)
    n0 = floor(u);  fr = u - n0;                                % neighbours n0, n0+1 and weights
    ok = n0 >= 0 & n0 + 1 <= Nt - 1;
    idx = find(ok);  m = numel(idx);
    Ip(k+1:k+m) = idx;        Jc(k+1:k+m) = (e-1)*Nt + n0(idx) + 1;  Vw(k+1:k+m) = we(e)*(1 - fr(idx));  k = k + m;
    Ip(k+1:k+m) = idx;        Jc(k+1:k+m) = (e-1)*Nt + n0(idx) + 2;  Vw(k+1:k+m) = we(e)*fr(idx);        k = k + m;
end
B = sparse(Ip(1:k), Jc(1:k), Vw(1:k), Np, Nt*Ne);
fprintf('B: %d x %d, nnz = %d (%.0f per row), density = %.3f %%\n', size(B,1), size(B,2), nnz(B), nnz(B)/Np, 100*nnz(B)/numel(B));
fprintf('memory: sparse ~%.1f MB vs dense %.1f GB (double)\n', nnz(B)*16/1e6, numel(B)*8/1e9);

%% 3. Simulate analytic channel data for three point targets (single scattering)
sig_f = frac_bw*f0/2/sqrt(2*log(2));  sig_t = 1/(2*pi*sig_f);    % Gaussian pulse, s
g = @(tt) exp(-tt.^2/(2*sig_t^2)) .* exp(2i*pi*f0*tt);          % analytic pulse
tg = [0 18; -1.5 20; 1.5 23]*1e-3;                              % targets [x z] (m)
d = zeros(Nt, Ne);
for s = 1:size(tg,1)
    for e = 1:Ne
        tau = (tg(s,2) + sqrt((tg(s,1) - xe(e))^2 + tg(s,2)^2))/c;
        d(:,e) = d(:,e) + g(t - tau);
    end
end
d = d + 0.01*(randn(Nt,Ne) + 1i*randn(Nt,Ne))/sqrt(2);         % weak white noise

%% 4. Beamform: one sparse matrix-vector product
img = B * d(:);                                 % complex image, Np x 1
I = reshape(abs(img), numel(zg), numel(xg));
IdB = 20*log10(I/max(I(:)));
% lateral -6 dB width at the first target
[~, iz] = min(abs(zg - 18e-3));  prof = I(iz,:)/max(I(iz,:));
[~, im] = max(prof);  l = im; while prof(l) > 0.5, l = l - 1; end;  r = im; while prof(r) > 0.5, r = r + 1; end
xl = xg(l) + (0.5 - prof(l))/(prof(l+1) - prof(l))*dxp;  xr = xg(r-1) + (0.5 - prof(r-1))/(prof(r) - prof(r-1))*dxp;
fprintf('lateral -6 dB width at z = 18 mm: %.2f mm (theory 1.21 lambda z / D = %.2f mm)\n', (xr - xl)*1e3, 1.21*lambda*18e-3/D*1e3);

%% 5. The adjoint B' maps an image back to channel data (back-projection)
dback = reshape(B' * img, Nt, Ne);              % every pixel re-emits its value along its own hyperbola
% dot test: <B d, g> must equal <d, B' g> for any d, g
dr = randn(Nt*Ne,1) + 1i*randn(Nt*Ne,1);  gr = randn(Np,1) + 1i*randn(Np,1);
lhs = gr' * (B*dr);  rhs = (B'*gr)' * dr;
fprintf('dot test |<Bd,g> - <d,B''g>| / |<Bd,g>| = %.1e\n', abs(lhs - rhs)/abs(lhs));

%% 6. B is H^H up to the pulse matched filter: compare B' e_j (pulse-filtered) with H e_j
[~, j] = min(abs(Xp(:) - 0) + abs(Zp(:) - 20e-3));            % pixel at (0, 20 mm)
aj = reshape(full(B(j,:)).', Nt, Ne);           % B' e_j: two-tap kernels on the hyperbola of pixel j
tp = (-ceil(4*sig_t*fs):ceil(4*sig_t*fs)).'/fs; gp = g(tp);
ajf = zeros(Nt, Ne);  hj = zeros(Nt, Ne);
for e = 1:Ne
    full_c = conv(aj(:,e), gp);  ajf(:,e) = full_c((1:Nt) + (numel(gp)-1)/2);   % 'same'-centred convolution
    tau = (Zp(j) + sqrt((Xp(j) - xe(e))^2 + Zp(j)^2))/c;  hj(:,e) = g(t - tau);  % H e_j: true echo of a point at pixel j
end
rho = abs(ajf(:)'*hj(:)) / (norm(ajf(:))*norm(hj(:)));
fprintf('|corr(pulse * B''e_j, H e_j)| = %.4f\n', rho);

%% 7. Figures
figure;
subplot(2,2,1); imagesc(1:Ne, t*1e6, real(d)); colorbar; xlabel('element'); ylabel('time (\mus)'); title('channel data d (real part)'); ylim([22 32]);
subplot(2,2,2); imagesc(xg*1e3, zg*1e3, IdB); caxis([-40 0]); axis image; colorbar;
xlabel('x (mm)'); ylabel('z (mm)'); title('image |B d| (dB)');
ix0 = find(abs(xg) < dxp/2);  rows0 = (ix0-1)*numel(zg) + (1:numel(zg));   % the 101 pixels of the x = 0 image column
subplot(2,2,3); spy(B(rows0, 1:4*Nt)); xlabel('column (e-1) N_t + n, elements 1-4'); ylabel('pixel (x = 0, z = 15..25 mm)');
title('spy(B): delay curves, 2 taps per element');
subplot(2,2,4); imagesc(1:Ne, t*1e6, abs(dback)); colorbar; xlabel('element'); ylabel('time (\mus)'); title('|B^H B d|: back-projection'); ylim([22 32]);
colormap(gray);
