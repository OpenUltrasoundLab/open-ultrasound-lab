%% m26_rayleigh_rician_k.m  --  envelope statistics: Rayleigh, Rician, K, and the log image
% Purpose : draw speckle IQ samples from three models and check their first-order
%           statistics against closed forms:
%             Rayleigh : circular complex Gaussian (fully developed speckle)
%             Rician   : coherent phasor A plus diffuse speckle of power s, r = A^2/s
%             K        : exponential speckle whose local mean fluctuates as Gamma(alpha)
%           For each: intensity contrast C_I = std(I)/mean(I), Nakagami m = 1/C_I^2,
%           envelope SNR0 = mean(R)/std(R), and the mean and std of the log image
%           10*log10(I/<I>) (Rayleigh: -2.51 dB and 5.57 dB).
% Change  : r_dB (coherent-to-diffuse ratios), alphas (K shape), Ns (samples).
% Expect  : Rayleigh C_I = 1, m = 1, SNR0 = 1.913, log mean -2.51 dB, log std 5.57 dB.
%           Rician  m = (1+r)^2/(1+2r): 1.33 (0 dB), 2.77 (6 dB), 5.76 (10 dB).
%           K       m = alpha/(alpha+2): 0.33 (alpha 1), 0.67 (alpha 4), 0.83 (alpha 10).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);
Ns = 400000;                                  % samples (independent resolution cells)
cg = @(n) (randn(1,n) + 1i*randn(1,n)) / sqrt(2);   % circular Gaussian, E|z|^2 = 1
stat = @(I) [std(I)/mean(I), mean(I)^2/var(I), mean(sqrt(I))/std(sqrt(I)), ...
             mean(10*log10(I/mean(I))), std(10*log10(I))];

%% 1. Rayleigh (fully developed speckle)
I_ray = abs(cg(Ns)).^2;
s = stat(I_ray);
fprintf('model            C_I     m      SNR0   logmean(dB) logstd(dB)   theory\n');
fprintf('Rayleigh       %6.3f %6.3f %6.3f   %7.2f   %7.2f      C_I=1 m=1 SNR0=%.3f mean=%.2f std=%.2f\n', ...
    s, sqrt(pi/(4-pi)), -10/log(10)*0.5772157, 10/log(10)*pi/sqrt(6));

%% 2. Rician: coherent component A on top of diffuse power s = 1
r_dB = [0 6 10];
I_ric = zeros(numel(r_dB), Ns);
for k = 1:numel(r_dB)
    r = 10^(r_dB(k)/10);
    z = sqrt(r) + cg(Ns);                     % A = sqrt(r), diffuse power 1
    I_ric(k,:) = abs(z).^2;
    s = stat(I_ric(k,:));
    fprintf('Rician r=%2d dB %6.3f %6.3f %6.3f   %7.2f   %7.2f      C_I=%.3f m=%.3f\n', ...
        r_dB(k), s, sqrt(1+2*r)/(1+r), (1+r)^2/(1+2*r));
end

%% 3. K distribution as a compound: I = x * E, x ~ Gamma(alpha, 1/alpha), E ~ Exp(1)
% Gamma variates without toolboxes: Marsaglia-Tsang for shape >= 1, boosted for shape < 1.
alphas = [1 4 10];
I_K = zeros(numel(alphas), Ns);
for k = 1:numel(alphas)
    a = alphas(k);  ab = a + (a < 1);         % boost shape below 1
    d = ab - 1/3;  c = 1/sqrt(9*d);
    x = zeros(1, Ns);  todo = 1:Ns;
    while ~isempty(todo)
        n = numel(todo);  g = randn(1, n);  v = (1 + c*g).^3;  u = rand(1, n);
        ok = v > 0 & log(u) < 0.5*g.^2 + d - d*v + d*log(max(v, eps));
        x(todo(ok)) = d * v(ok);  todo = todo(~ok);
    end
    if a < 1, x = x .* rand(1, Ns).^(1/a); end
    x = x / a;                                % mean 1, variance 1/alpha
    I_K(k,:) = x .* (-log(rand(1, Ns)));      % exponential speckle with fluctuating mean
    s = stat(I_K(k,:));
    fprintf('K alpha=%-4g   %6.3f %6.3f %6.3f   %7.2f   %7.2f      C_I=%.3f m=%.3f\n', ...
        a, s, sqrt(1+2/a), a/(a+2));
end

%% 4. Figures: envelope histograms with closed-form pdfs; log-image histograms
Rgrid = linspace(0.005, 4, 400);
figure('Position', [80 80 1150 380]);
subplot(1,3,1);                               % Rayleigh and Rician envelopes
cols = [0.2 0.4 0.8; 0.8 0.3 0.5; 0.7 0.5 0.1; 0.1 0.55 0.4];
edges = linspace(0, 6, 61);  ctr = (edges(1:end-1) + edges(2:end))/2;  bw = edges(2) - edges(1);
R = sqrt(I_ray);  cnt = histc(R, edges);  plot(ctr, cnt(1:end-1)/(Ns*bw), 'o', 'Color', cols(1,:)); hold on;
plot(Rgrid, 2*Rgrid.*exp(-Rgrid.^2), '-', 'Color', cols(1,:), 'LineWidth', 1.5);
Rg2 = linspace(0.005, 6, 500);
for k = 1:numel(r_dB)
    r = 10^(r_dB(k)/10);  A = sqrt(r);
    R = sqrt(I_ric(k,:));  cnt = histc(R, edges);
    plot(ctr, cnt(1:end-1)/(Ns*bw), 'o', 'Color', cols(k+1,:));
    % Rice pdf with diffuse power 1: 2R exp(-(R^2+A^2)) I0(2AR), via scaled besseli
    plot(Rg2, 2*Rg2.*exp(-(Rg2 - A).^2).*besseli(0, 2*A*Rg2, 1), '-', 'Color', cols(k+1,:), 'LineWidth', 1.5);
end
xlabel('envelope R (diffuse rms = 1)'); ylabel('pdf');
title('Rayleigh and Rician (r = 0, 6, 10 dB)'); xlim([0 6]);

subplot(1,3,2);                               % K envelopes (mean intensity 1) with Nakagami fit
edges = linspace(0, 4, 81);  ctr = (edges(1:end-1) + edges(2:end))/2;  bw = edges(2) - edges(1);
for k = 1:numel(alphas)
    a = alphas(k);  b = 2*sqrt(a);            % <I> = 1
    R = sqrt(I_K(k,:));  cnt = histc(R, edges);
    semilogy(ctr, cnt(1:end-1)/(Ns*bw), 'o', 'Color', cols(k+1,:)); hold on;
    pK = 2*b/gamma(a) * (b*Rgrid/2).^a .* besselk(a-1, b*Rgrid);
    semilogy(Rgrid, pK, '-', 'Color', cols(k+1,:), 'LineWidth', 1.5);
    m = a/(a+2);                              % Nakagami with the same m and Omega = 1
    semilogy(Rgrid, 2*m^m/gamma(m) * Rgrid.^(2*m-1) .* exp(-m*Rgrid.^2), '--', 'Color', cols(k+1,:));
end
semilogy(Rgrid, 2*Rgrid.*exp(-Rgrid.^2), 'k-', 'LineWidth', 1);
ylim([1e-4 3]); xlabel('envelope R (<R^2> = 1)'); ylabel('pdf (log scale)');
title('K: alpha = 1, 4, 10 (solid), Nakagami (dashed)');

subplot(1,3,3);                               % log-compressed speckle
edges = -45:1:15;  ctr = (edges(1:end-1) + edges(2:end))/2;
D = 10*log10(I_ray);  cnt = histc(D, edges);  plot(ctr, cnt(1:end-1)/Ns, 'o', 'Color', cols(1,:)); hold on;
Dg = linspace(-45, 15, 400);  Ig = 10.^(Dg/10);
plot(Dg, log(10)/10 * Ig .* exp(-Ig), '-', 'Color', cols(1,:), 'LineWidth', 1.5);
D = 10*log10(I_ric(3,:)/mean(I_ric(3,:)));  cnt = histc(D, edges);  plot(ctr, cnt(1:end-1)/Ns, 's', 'Color', cols(4,:));
D = 10*log10(I_K(1,:));  cnt = histc(D, edges);  plot(ctr, cnt(1:end-1)/Ns, '^', 'Color', cols(2,:));
xlabel('10 log_{10}(I/<I>)  (dB)'); ylabel('fraction per dB');
title('log image: Rayleigh (line), Rician 10 dB, K alpha=1'); xlim([-45 15]);
