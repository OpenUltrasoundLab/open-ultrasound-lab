%% m29_ultrafast_power_doppler.m
% Purpose : Ultrafast power Doppler on a simulated slow-time ensemble.
%           Clutter: speckle tissue that moves axially (3.3 Hz, 50 um, with a
%           depth gradient so the motion is not rigid) at 40 dB above blood.
%           Blood: five small vessels (150 um wide) with axial speeds from
%           1 to 30 mm/s, each pixel a complex Gaussian process with a
%           Gaussian Doppler spectrum. White noise 6 dB above blood.
%           Part 1: SVD clutter filter on the Casorati matrix (pixels x time),
%                   clutter rank Kc = singular powers > 10 dB above the median.
%           Part 2: power Doppler CNR per vessel versus ensemble length N;
%                   theory: CNR = (P_b / P_n) sqrt(N - Kc) for independent noise.
%           Part 3: SVD versus an ideal FFT high-pass (|f| > 75 Hz) at N = 200.
% Change  : f0, PRF, motion amplitude Am, CBRdB, BNRdB, fcut, the vessel speeds.
% Expect  : f0 = 15 MHz, PRF = 1 kHz, peak tissue Doppler 20 Hz (30 Hz at the
%           bottom). Kc = 4 at N = 400; the first 4 components hold 99.9993 %
%           of the tissue energy. CNR follows sqrt(N - Kc): 2.47 -> 5.00 from
%           N = 100 to 400 (measured ratio 1.98). Slow vessels fall below theory
%           at N = 25-50 (blood not yet decorrelated: it shares the tissue
%           subspace). The FFT high-pass leaks tissue through the rectangular
%           window (motion not periodic in the ensemble): background 11.5 dB
%           above noise, no vessel detectable, while SVD leaves 0.99 x noise.
rng(1);

%% Parameters
c = 1540; f0 = 15e6; PRF = 1e3; T = 1/PRF;     % compounded frame rate (Hz)
Nmax = 400;                                     % longest ensemble (frames)
nz = 64; nx = 64; dpix = 50e-6;                 % pixels and pixel size (m)
P = nz*nx;  t = (0:Nmax-1)*T;
CBRdB = 40;  BNRdB = -6;                        % clutter/blood and blood/noise per sample (dB)
Am = 50e-6;  fm = 3.3;  grad = 0.5;             % motion amplitude (m), rate (Hz, not periodic in the window), depth gradient
fcut = 75;                                      % FFT high-pass cutoff for comparison (Hz)
vax = [1 3 10 20 30]*1e-3;                      % axial blood speeds (m/s)
zc  = [8 20 32 44 56];                          % vessel centre rows (pixels)
fD  = 2*vax*f0/c;
fprintf('lambda = %.1f um; Doppler of the five vessels: %s Hz (Nyquist %d Hz)\n', c/f0*1e6, mat2str(round(fD*10)/10), PRF/2);
fprintf('Peak tissue velocity %.2f mm/s -> peak tissue Doppler %.1f Hz (x%.1f at the bottom)\n', 2*pi*fm*Am*1e3, 2*2*pi*fm*Am*f0/c, 1 + grad);

%% Signals
[XX, ZZ] = meshgrid(1:nx, 1:nz);
lab = zeros(nz, nx);                             % vessel label per pixel (0 = background)
for k = 1:5, lab(abs(ZZ - zc(k) - 0.15*(XX - nx/2)) <= 1.5) = k; end   % slightly oblique, 3 px wide
% Tissue: fixed speckle amplitude x motion phase 4 pi f0 d(z,t)/c
depthw = 1 + grad*(ZZ(:) - 1)/(nz - 1);          % P x 1 displacement weight
d = Am*sin(2*pi*fm*t + 0.7);                          % 1 x N displacement (m)
spk = (randn(P,1) + 1j*randn(P,1))/sqrt(2);
tiss = spk .* exp(1j*4*pi*f0/c * depthw * d);    % P x N, unit power
% Blood: white noise shaped by a Gaussian spectrum centred at fD (per vessel)
Lf = 1024;  ff = ((0:Lf-1) - Lf/2)/Lf*PRF;  blood = zeros(P, Nmax);
for k = 1:5
    idx = find(lab(:) == k);  sb = 20 + 0.3*fD(k);                % spectral std (Hz)
    fa = mod(ff - fD(k) + PRF/2, PRF) - PRF/2;                     % aliased offset from fD
    H = ifftshift(exp(-fa.^2/(4*sb^2)));                           % amplitude filter -> power std sb
    W = (randn(numel(idx), Lf) + 1j*randn(numel(idx), Lf))/sqrt(2);
    B = ifft(bsxfun(@times, fft(W, [], 2), H), [], 2);
    B = B(:, 1:Nmax);  blood(idx, :) = B / sqrt(mean(abs(B(:)).^2));
end
noise = (randn(P, Nmax) + 1j*randn(P, Nmax))/sqrt(2) * 10^(-BNRdB/20);
S = 10^(CBRdB/20)*tiss + blood + noise;          % Casorati matrix, P x N

%% Part 1: SVD of the full ensemble
[U, Sg, V] = svd(S, 'econ');  sv = diag(Sg);
s2 = sv.^2;  floorP = median(s2);
Kc = find(s2 > 10*floorP, 1, 'last');            % clutter: > 10 dB above the median singular power
fprintf('Part 1: N = %d, clutter rank Kc = %d (singular power > 10 dB above median); s1/s_median = %.1f dB\n', Nmax, Kc, 10*log10(s2(1)/floorP));
fprintf('        tissue energy in first Kc components: %.4f %%  of total tissue energy left = %.2e\n', ...
    100*(1 - norm(10^(CBRdB/20)*tiss - U(:,1:Kc)*(U(:,1:Kc)'*10^(CBRdB/20)*tiss), 'fro')^2/norm(10^(CBRdB/20)*tiss, 'fro')^2), ...
    norm(10^(CBRdB/20)*tiss - U(:,1:Kc)*(U(:,1:Kc)'*10^(CBRdB/20)*tiss), 'fro')^2/norm(10^(CBRdB/20)*tiss, 'fro')^2);

%% Part 2: power Doppler CNR vs ensemble length
Ns = [25 50 100 200 400];  bg = lab(:) == 0;  cnr = zeros(numel(Ns), 5);  Kn = zeros(size(Ns));
fprintf('Part 2: CNR = (mean PD vessel - mean PD background) / std PD background\n');
fprintf('   N    Kc   V1(1mm/s) V2(3)  V3(10)  V4(20)  V5(30)   theory (Pb/Pn) sqrt(N - Kc)\n');
for i = 1:numel(Ns)
    N = Ns(i);  Sn = S(:, 1:N);
    [Un, Sgn, Vn] = svd(Sn, 'econ');  q = diag(Sgn).^2;  K = find(q > 10*median(q), 1, 'last');  Kn(i) = K;
    F = Sn - Un(:,1:K)*(Sgn(1:K,1:K)*Vn(:,1:K)');
    PD = mean(abs(F).^2, 2);
    for k = 1:5, cnr(i,k) = (mean(PD(lab(:)==k)) - mean(PD(bg))) / std(PD(bg)); end
    fprintf('  %3d   %2d   %5.2f   %5.2f   %5.2f   %5.2f   %5.2f    %5.2f\n', N, K, cnr(i,:), 10^(BNRdB/10)*sqrt(N - K));
    if N == 200, PDsvd = PD; F200 = F; end
end
fprintf('CNR ratio N = 400 vs 100 (mean of vessels): %.2f (sqrt(4) = 2)\n', mean(cnr(5,:)./cnr(3,:)));

%% Part 3: SVD vs FFT high-pass at N = 200
N = 200;  Sn = S(:, 1:N);
fk = ((0:N-1) - N*(0:N-1 >= N/2))*PRF/N;         % FFT bin frequencies (Hz)
G = fft(Sn, [], 2);  G(:, abs(fk) <= fcut) = 0;  Fh = ifft(G, [], 2);
PDhp = mean(abs(Fh).^2, 2);
cnrH = zeros(1,5);  cnrS = zeros(1,5);
for k = 1:5
    cnrH(k) = (mean(PDhp(lab(:)==k)) - mean(PDhp(bg))) / std(PDhp(bg));
    cnrS(k) = (mean(PDsvd(lab(:)==k)) - mean(PDsvd(bg))) / std(PDsvd(bg));
end
resH = mean(PDhp(bg)) / (10^(-BNRdB/10));   resS = mean(PDsvd(bg)) / (10^(-BNRdB/10));
fprintf('Part 3 (N = 200): background PD / noise power: SVD %.2f, high-pass %.1f (%.1f dB above noise)\n', resS, resH, 10*log10(resH));
fprintf('   vessel:   1 mm/s  3 mm/s  10 mm/s  20 mm/s  30 mm/s\n   SVD CNR  %6.2f  %6.2f  %6.2f   %6.2f   %6.2f\n   HP  CNR  %6.2f  %6.2f  %6.2f   %6.2f   %6.2f\n', cnrS, cnrH);

%% Figures
figure('Position', [40 40 1150 640]);
subplot(2,3,1); plot(10*log10(s2/s2(1)), '.-'); hold on; plot([Kc Kc], [-80 0], 'k--'); grid on;
xlabel('singular index'); ylabel('s_k^2 / s_1^2 (dB)'); title(sprintf('Singular spectrum, N = %d, Kc = %d', Nmax, Kc)); xlim([1 Nmax]);
subplot(2,3,2); plot(t(1:200)*1e3, real(V(1:200,1)), t(1:200)*1e3, real(V(1:200,Kc+5)), 'LineWidth', 1.2); grid on;
xlabel('slow time (ms)'); ylabel('Re v_k'); legend('k = 1 (tissue)', sprintf('k = %d (blood/noise)', Kc+5)); title('Temporal singular vectors');
subplot(2,3,3); loglog(Ns, cnr, 'o-'); hold on; loglog(Ns, 10^(BNRdB/10)*sqrt(Ns), 'k--', 'LineWidth', 1.5); grid on;
xlabel('ensemble length N'); ylabel('power Doppler CNR'); title('Sensitivity vs ensemble length'); legend('1', '3', '10', '20', '30 mm/s', 'sqrt(N) P_b/P_n', 'Location', 'northwest');
subplot(2,3,4); imagesc((1:nx)*dpix*1e3, (1:nz)*dpix*1e3, reshape(10*log10(mean(abs(S(:,1:200)).^2,2)), nz, nx)); axis image; colorbar; colormap(gray);
title('Unfiltered mean power (dB)'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,5); imagesc((1:nx)*dpix*1e3, (1:nz)*dpix*1e3, reshape(10*log10(PDsvd*10^(BNRdB/10)), nz, nx), [-2 3]); axis image; colorbar;
title('Power Doppler, SVD, N = 200 (dB re noise)'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,6); imagesc((1:nx)*dpix*1e3, (1:nz)*dpix*1e3, reshape(10*log10(PDhp*10^(BNRdB/10)), nz, nx), [-2 20]); axis image; colorbar;
title(sprintf('PD, high-pass > %d Hz (dB re noise)', fcut)); xlabel('x (mm)'); ylabel('z (mm)');
