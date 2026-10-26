%% m07_complex_gaussian_speckle.m  --  speckle IQ is circular complex Gaussian
% Purpose : build IQ samples as sums of many random phasors (sub-resolution scatterers),
%           then check every property of a circular complex Gaussian: I and Q
%           uncorrelated with equal variance, zero pseudo-covariance E[z^2],
%           Rayleigh envelope (mean/std = 1.91), exponential intensity (mean/std = 1),
%           uniform phase. Then break circularity with I/Q imbalance, and show that
%           I and |z| are uncorrelated but dependent.
% Change  : Nscat (scatterers per cell), gain_err, skew_deg, A_coh (coherent component).
% Expect  : |P|/R ~ 0.01 (1/sqrt(Ns) floor), envelope SNR 1.91, intensity SNR 1.00,
%           corr(I,|z|) ~ 0, corr(I^2,|z|^2) ~ 0.707; with imbalance |P|/R ~ 0.42
%           while the envelope SNR barely moves.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);

%% 1. Random phasor sums (central limit theorem)
Ns    = 20000;                  % number of IQ samples (independent resolution cells)
Nscat = 30;                     % scatterers per resolution cell
amp   = rand(Nscat, Ns);                          % arbitrary positive amplitudes
phs   = 2*pi*rand(Nscat, Ns);                     % independent uniform phases
z     = sum(amp .* exp(1i*phs), 1) / sqrt(Nscat); % 1 x Ns complex samples
z     = z / sqrt(mean(abs(z).^2));                % normalise E|z|^2 = 1
I = real(z);  Q = imag(z);

%% 2. Second-order statistics
C_IQ = [mean(I.^2) mean(I.*Q); mean(I.*Q) mean(Q.^2)];   % real 2x2 covariance (zero mean)
R    = mean(abs(z).^2);                                  % E[z z*]
P    = mean(z.^2);                                       % pseudo-covariance E[z z]
fprintf('real covariance of (I,Q) = [%.3f %.3f; %.3f %.3f]\n', C_IQ');
fprintf('R = E|z|^2 = %.3f   P = E[z^2] = %.4f %+.4fj   |P|/R = %.4f (floor 1/sqrt(Ns) = %.4f)\n', ...
    R, real(P), imag(P), abs(P)/R, 1/sqrt(Ns));

%% 3. Envelope, intensity, phase
a  = abs(z);  Iint = a.^2;  ph = angle(z);
fprintf('envelope mean/std = %.3f (Rayleigh sqrt(pi/(4-pi)) = %.3f)\n', mean(a)/std(a), sqrt(pi/(4-pi)));
fprintf('intensity mean/std = %.3f (exponential: 1)\n', mean(Iint)/std(Iint));
fprintf('E|z| = %.3f (theory sqrt(pi)/2 = %.3f)\n', mean(a), sqrt(pi)/2);

%% 4. Uncorrelated is not independent
cc = @(u, v) mean((u-mean(u)).*(v-mean(v))) / (std(u,1)*std(v,1));
fprintf('corr(I, Q) = %.4f   corr(I, |z|) = %.4f   corr(I^2, |z|^2) = %.3f (theory 1/sqrt(2) = 0.707)\n', ...
    cc(I, Q), cc(I, a), cc(I.^2, Iint));

%% 5. Break circularity: I/Q gain and phase imbalance in the demodulator
gain_err = 0.2;  skew_deg = 10;  A_coh = 0;       % fractional gain error, phase skew (deg), coherent mean
Ib = (1 + gain_err) * I + A_coh;
Qb = (1 - gain_err) * (cos(skew_deg*pi/180) * Q + sin(skew_deg*pi/180) * I);
zb = Ib + 1i*Qb;  zb0 = zb - mean(zb);
Pb = mean(zb0.^2);  Rb = mean(abs(zb0).^2);
fprintf('with imbalance: |P|/R = %.3f   envelope mean/std = %.3f   phase histogram max/min = %.2f\n', ...
    abs(Pb)/Rb, mean(abs(zb))/std(abs(zb)), max(hist(angle(zb), 36)) / min(hist(angle(zb), 36)));

%% 6. Figures
figure('Position', [100 100 1100 360]);
subplot(1,3,1);
plot(I(1:3000), Q(1:3000), '.', 'MarkerSize', 3); hold on;
plot(Ib(1:3000), Qb(1:3000), '.', 'MarkerSize', 3, 'Color', [0.85 0.4 0.2]);
axis equal; axis([-3 3 -3 3]); grid on;
xlabel('I (a.u.)'); ylabel('Q (a.u.)'); title('IQ plane: circular (blue) vs imbalanced (orange)');
subplot(1,3,2);
centers = linspace(0.05, 3.45, 35);  bw = centers(2) - centers(1);
cnt = hist(a, centers);
bar(centers, cnt / (Ns*bw), 1, 'FaceColor', [0.6 0.7 0.9], 'EdgeColor', 'none'); hold on;
aa = linspace(0, 3.5, 300);
plot(aa, 2*aa/R .* exp(-aa.^2/R), 'r-', 'LineWidth', 1.8);
xlabel('envelope |z|'); ylabel('density'); title('envelope vs Rayleigh');
subplot(1,3,3);
centers = linspace(-pi + pi/36, pi - pi/36, 36);  bw = centers(2) - centers(1);
cnt = hist(ph, centers);
bar(centers, cnt / (Ns*bw), 1, 'FaceColor', [0.6 0.7 0.9], 'EdgeColor', 'none'); hold on;
plot([-pi pi], [1 1]/(2*pi), 'r-', 'LineWidth', 1.8);
xlim([-pi pi]); ylim([0 0.25]);
xlabel('phase arg z (rad)'); ylabel('density'); title('phase vs uniform 1/(2\pi)');
