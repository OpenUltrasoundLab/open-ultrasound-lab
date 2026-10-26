%% m06_tsvd_picard.m  --  Module 06: truncated SVD, Tikhonov and the Picard plot
% Purpose : Deconvolve a blurred, noisy 1-D reflectivity profile b = A x + noise,
%           where A is a 64 x 64 Gaussian blur (the same problem as Demo 06.5).
%           Compare the plain pseudo-inverse, truncated SVD (TSVD) and Tikhonov,
%           read the Picard plot (sigma_i, |u_i' b|, |u_i' b|/sigma_i) and see
%           why the truncation index sits where |u_i' b| meets the noise floor.
% Change  : blur width sg (samples), noise levels nl (dB re signal RMS), n.
% Expect  : kappa(A) = 5.7e12 (sigma_64/sigma_1 = 1.75e-13); MATLAB's default
%           pinv tolerance (7e-15) keeps all 64 terms, so pinv(A)*b has relative
%           error ~1e10 at -40 dB noise. Best TSVD k = 27 (error 0.40), and the
%           Picard rule "last i with |u_i' b0| above the noise std" predicts 26;
%           Tikhonov at its best alpha gives 0.40 too. At -60 dB the best k moves
%           only to 31: sigma_i of a Gaussian blur falls super-exponentially.
% Runtime : < 1 s.
clear; close all; rng(1);

%% 1. Forward model
n = 64; sg = 2.5;                              % samples, Gaussian blur std (samples)
[I, J] = ndgrid(1:n, 1:n);
A = exp(-(I - J).^2/(2*sg^2))/(sqrt(2*pi)*sg);  % symmetric Toeplitz blur, rows sum to ~1
x = zeros(n,1);
x(11:18) = 0.6;                                % a uniform layer
x(39:50) = 0.35 + 0.02*((39:50) - 39);         % a ramped layer
x(28) = 1; x(31) = -0.7;                       % two point reflectors 3 samples apart
b0 = A*x;

[U, S, V] = svd(A); s = diag(S);
fprintf('kappa(A) = sigma_1/sigma_n = %.2e  (sigma_n/sigma_1 = %.2e)\n', s(1)/s(end), s(end)/s(1));
fprintf('default rank/pinv tolerance max(m,n)*eps*sigma_1 = %.1e -> numerical rank %d of %d\n', ...
    max(size(A))*eps(s(1)), sum(s > max(size(A))*eps(s(1))), n);

%% 2. Spectral filters at two noise levels
nl = [-40 -60];                                % noise re signal RMS (dB)
rmsb = sqrt(mean(b0.^2));
fprintf('\n noise   pinv err   best TSVD k  err     best Tikhonov alpha  err\n');
for q = 1:numel(nl)
    sn = rmsb*10^(nl(q)/20);                   % noise std per sample
    b = b0 + sn*randn(n,1);
    beta = U'*b;                               % data coefficients u_i' b
    errK = zeros(n,1);
    for k = 1:n
        xk = V(:,1:k)*(beta(1:k)./s(1:k));
        errK(k) = norm(xk - x)/norm(x);
    end
    [eK, kB] = min(errK);
    alphas = logspace(-6, 0, 121); errA = zeros(size(alphas));
    for j = 1:numel(alphas)
        f = s.^2./(s.^2 + alphas(j)^2);        % Tikhonov filter factors
        errA(j) = norm(V*(f.*beta./s) - x)/norm(x);
    end
    [eA, jA] = min(errA);
    ep = norm(pinv(A)*b - x)/norm(x);
    fprintf(' %3d dB   %8.1e      %2d     %.3f        %.1e          %.3f\n', nl(q), ep, kB, eK, alphas(jA), eA);
    % Picard prediction: keep terms while |u_i' b0| > noise std
    kP = find(abs(U'*b0) > sn, 1, 'last');
    fprintf('          last i with |u_i''b0| > noise std: %d\n', kP);
    if q == 1, b1 = b; beta1 = beta; errK1 = errK; kB1 = kB; sn1 = sn; end
end

%% 3. Figures (first noise level)
figure('Position', [60 60 1000 380]);
subplot(1,2,1);
semilogy(1:n, s, 'k.-', 1:n, abs(beta1), 'mo', 1:n, abs(beta1)./s, 'b.'); hold on;
semilogy([1 n], sn1*[1 1], 'm--'); semilogy([kB1 kB1]+0.5, [1e-16 1e8], 'r-');
ylim([1e-16 1e8]); grid on; xlabel('index i'); ylabel('magnitude');
legend('\sigma_i', '|u_i^T b|', '|u_i^T b| / \sigma_i', 'noise std', 'best k', 'Location', 'southwest');
title(sprintf('Picard plot, noise %d dB', nl(1)));
subplot(1,2,2);
xk = V(:,1:kB1)*(beta1(1:kB1)./s(1:kB1));
plot(1:n, x, 'k--', 1:n, b1, 'm', 1:n, xk, 'b', 'LineWidth', 1.4); grid on;
xlabel('depth sample'); ylabel('reflectivity');
legend('true x', 'data b', sprintf('TSVD, k = %d', kB1), 'Location', 'northeast');
title(sprintf('TSVD estimate, relative error %.2f', errK1(kB1)));
