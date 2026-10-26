%% m12_crlb_amplitude_mle.m
% Purpose : the Cramer-Rao bound in its simplest settings.
%   Part 1: DC level A in N noise samples. Sample mean, median and midrange
%           under Gaussian, Laplace and uniform noise of equal variance.
%           Prints N*var/sigma^2, so the Gaussian CRLB is 1 and the Laplace
%           CRLB is 0.5. Uniform noise violates the regularity condition
%           (its support moves with A): no CRLB, and the midrange beats 1/N.
%   Part 2: amplitude of a known pulse s in complex circular Gaussian noise
%           (IQ data). Real amplitude: var >= sigma^2/(2||s||^2).
%           Complex amplitude: E|a_hat - a|^2 >= sigma^2/||s||^2.
%           Same pulse in real noise of variance sigma^2: sigma^2/||s||^2.
%   Part 3: MLE invariance. The MLE of |a|^2 is |a_hat|^2, which is biased
%           upward by sigma^2/||s||^2.
% Change  : N, nTrial. Expect Laplace-median N*var to approach 0.5 slowly
%           (0.58 at N = 101).
% Expect  : N = 25: Gaussian mean 1.00, median ~1.54, midrange ~3.3;
%           Laplace mean 1.00, median ~0.67; uniform midrange ~0.214
%           (theory 2*3*N/((N+1)(N+2)) = 0.2137). Part 2 ratios ~1.00.
rng(1);

%% Part 1: DC level, three noise laws with unit variance
N = 25;            % samples per estimate
nTrial = 20000;    % Monte Carlo trials
A = 0;             % true level (any value; estimators are shift-equivariant)
lawName = {'Gaussian', 'Laplace ', 'uniform '};
crlbN = [1, 0.5, NaN];   % N*CRLB/sigma^2 (none for uniform)
res = zeros(3, 3);       % rows: law, cols: mean, median, midrange
estAll = cell(3, 1);
for law = 1:3
    switch law
        case 1
            w = randn(nTrial, N);
        case 2   % Laplace with variance 1: scale b = 1/sqrt(2)
            u = rand(nTrial, N) - 0.5;
            w = -(1/sqrt(2)) * sign(u) .* log(1 - 2*abs(u));
        case 3   % uniform on [-sqrt(3), sqrt(3)]: variance 1
            w = (rand(nTrial, N) - 0.5) * 2*sqrt(3);
    end
    y = A + w;
    est = [mean(y, 2), median(y, 2), (max(y, [], 2) + min(y, [], 2))/2];
    estAll{law} = est;
    res(law, :) = N * var(est - A);
end
fprintf('Part 1: N*var/sigma^2 for N = %d (%d trials)\n', N, nTrial);
fprintf('  noise      mean    median  midrange   N*CRLB\n');
for law = 1:3
    fprintf('  %s  %6.3f  %6.3f  %7.3f   %5.2f\n', lawName{law}, res(law, :), crlbN(law));
end
fprintf('  theory: median/Gaussian -> pi/2 = %.3f; Laplace median -> 0.5; uniform midrange = %.4f\n', ...
    pi/2, 6*N/((N+1)*(N+2)));

% N*var versus N for the Gaussian and Laplace cases
Ns = [3 5 9 17 33 65 129];
effN = zeros(numel(Ns), 4);  % Gauss mean, Gauss median, Laplace mean, Laplace median
for k = 1:numel(Ns)
    n = Ns(k); g = randn(5000, n);
    u = rand(5000, n) - 0.5; l = -(1/sqrt(2)) * sign(u) .* log(1 - 2*abs(u));
    effN(k, :) = n * [var(mean(g, 2)), var(median(g, 2)), var(mean(l, 2)), var(median(l, 2))];
end

%% Part 2: amplitude of a known complex pulse (IQ), real vs complex noise
fs = 40e6; f0 = 5e6; t = (-40:40)'/fs;          % time axis (s)
s = exp(-t.^2/(2*(0.15e-6)^2)) .* exp(1j*2*pi*f0*t);  % known complex pulse
sigma2 = 0.5;                                   % noise variance E|w|^2
a = 0.8;                                        % true real amplitude
M = 20000; E = sum(abs(s).^2);
W = sqrt(sigma2/2) * (randn(numel(s), M) + 1j*randn(numel(s), M));
Y = a*s + W;
proj = (s' * Y) / E;                            % s^H y / ||s||^2 for each trial
aReal = real(proj);                             % MLE when a is real (phase known)
aCplx = proj;                                   % MLE when a is complex
fprintf('\nPart 2: known pulse, ||s||^2 = %.2f, sigma^2 = %.2f\n', E, sigma2);
fprintf('  real a, complex noise : var = %.3e, CRLB sigma^2/(2||s||^2) = %.3e, ratio %.3f\n', ...
    var(aReal), sigma2/(2*E), var(aReal)/(sigma2/(2*E)));
fprintf('  complex a             : E|err|^2 = %.3e, CRLB sigma^2/||s||^2 = %.3e, ratio %.3f\n', ...
    mean(abs(aCplx - a).^2), sigma2/E, mean(abs(aCplx - a).^2)/(sigma2/E));
sr = real(s); Er = sum(sr.^2);                  % real pulse in real noise of variance sigma^2
yr = a*sr + sqrt(sigma2)*randn(numel(sr), M);
aR = (sr' * yr) / Er;
fprintf('  real pulse, real noise: var = %.3e, CRLB sigma^2/||s||^2 = %.3e, ratio %.3f\n', ...
    var(aR), sigma2/Er, var(aR)/(sigma2/Er));

%% Part 3: invariance is not unbiasedness
p = abs(aCplx).^2;                              % MLE of the power |a|^2
fprintf('\nPart 3: E[|a_hat|^2] = %.4f, |a|^2 = %.4f, bias = %.4f, predicted sigma^2/||s||^2 = %.4f\n', ...
    mean(p), a^2, mean(p) - a^2, sigma2/E);

%% Figures
figure('Position', [80 80 1000 380]);
subplot(1, 2, 1);
est = estAll{1}; edges = linspace(-1, 1, 61); ctr = (edges(1:end-1) + edges(2:end))/2;
h1 = histc(est(:, 1), edges); h2 = histc(est(:, 2), edges); h3 = histc(est(:, 3), edges);
bw = edges(2) - edges(1);
plot(ctr, h1(1:end-1)/(nTrial*bw), 'b-', ctr, h2(1:end-1)/(nTrial*bw), 'r-', ...
     ctr, h3(1:end-1)/(nTrial*bw), 'Color', [0.8 0.5 0], 'LineWidth', 1.4); hold on;
plot(ctr, sqrt(N/(2*pi)) * exp(-N*ctr.^2/2), 'k--', 'LineWidth', 1.2);
xlabel('estimate - A (noise std units)'); ylabel('density');
title(sprintf('Gaussian noise, N = %d', N));
legend('mean', 'median', 'midrange', 'CRLB pdf N(0,\sigma^2/N)');
subplot(1, 2, 2);
semilogx(Ns, effN(:, 1), 'b-o', Ns, effN(:, 2), 'r-o', Ns, effN(:, 3), 'b--s', Ns, effN(:, 4), 'r--s'); hold on;
semilogx(Ns, ones(size(Ns)), 'k:', Ns, 0.5*ones(size(Ns)), 'k-.');
xlabel('N (samples)'); ylabel('N \cdot var / \sigma^2'); ylim([0 2]);
title('Efficiency: Gaussian (solid) and Laplace (dashed) noise');
legend('mean, Gauss', 'median, Gauss', 'mean, Laplace', 'median, Laplace', 'CRLB Gauss', 'CRLB Laplace');
