%% m26_nakagami_density.m  --  how many scatterers make speckle "fully developed"?
% Purpose : sum unit phasors with uniform random phases in each resolution cell and
%           measure the envelope statistics as the number of scatterers per cell grows.
%             Poisson count, mean N  : C_I^2 = 1 + 1/N  ->  Nakagami m = N/(N+1)  (K-like, m < 1)
%             exactly N per cell     : C_I^2 = 1 - 1/N  ->  m = N/(N-1)          (sub-Rayleigh, m > 1)
%           The Poisson result also holds for PSF-weighted scatterers if N = n*Vc, where
%           n is the number density and Vc = (int |h|^2)^2 / int |h|^4 the correlation area.
%           Part 2 checks the spread of the moment estimator m^ = (E[R^2])^2 / Var(R^2)
%           from n independent cells (std about 2/sqrt(n) at m = 1).
% Change  : Nlist (mean scatterers per cell), Ncell, nlist (cells per estimate).
% Expect  : Poisson N = 10 -> m = 0.909, SNR0 = 1.86 (97 % of 1.913);
%           fixed N = 10 -> m = 1.111; std(m^) * sqrt(n) close to 2.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);
Ncell = 40000;                                 % independent resolution cells per point
Nlist = [1 2 3 5 10 20 50 100];                % mean number of scatterers per cell
poissrnd0 = @(lam, n) sum(cumsum(-log(rand(ceil(lam + 7*sqrt(lam) + 12), n)) / lam, 1) <= 1, 1); % Poisson via exponential gaps
nak  = @(R) mean(R.^2)^2 / var(R.^2);          % Nakagami m, moment estimator
snr0 = @(R) mean(R) / std(R);                  % envelope SNR (Rayleigh: 1.913)

%% 1. Statistics vs mean scatterer count
mP = zeros(size(Nlist));  sP = mP;  mF = mP;  sF = mP;
fprintf('  N    Poisson: m  (theory)  SNR0  | fixed N: m  (theory)  SNR0\n');
for k = 1:numel(Nlist)
    N = Nlist(k);
    cnt = poissrnd0(N, Ncell);                % Poisson number of scatterers in each cell
    maxc = max(cnt);
    ph = 2*pi*rand(maxc, Ncell);
    use = bsxfun(@le, (1:maxc)', cnt);        % keep the first cnt(j) phasors of column j
    RP = abs(sum(exp(1i*ph) .* use, 1));
    RF = abs(sum(exp(1i*2*pi*rand(N, Ncell)), 1));   % exactly N phasors per cell
    mP(k) = nak(RP);  sP(k) = snr0(RP);  mF(k) = nak(RF);  sF(k) = snr0(RF);
    if N == 1                                 % one unit phasor: R = 1 exactly, Var(R^2) = 0
        fprintf('%4d   %8.3f  (%5.3f)  %5.3f  |   R = 1 exactly (no speckle)\n', N, mP(k), N/(N+1), sP(k));
    else
        fprintf('%4d   %8.3f  (%5.3f)  %5.3f  | %8.3f  (%5.3f)  %5.3f\n', N, mP(k), N/(N+1), sP(k), mF(k), N/(N-1), sF(k));
    end
end
fprintf('Rayleigh limit: m = 1, SNR0 = %.3f\n', sqrt(pi/(4-pi)));

%% 2. Spread of the Nakagami estimate from n independent Rayleigh cells
nlist = [25 100 400 1600];  Ntrial = 4000;
fprintf('\n   n   mean(m^)   std(m^)   std*sqrt(n)\n');
for n = nlist
    R = abs(randn(n, Ntrial) + 1i*randn(n, Ntrial));     % Rayleigh envelopes
    I = R.^2;
    mh = mean(I, 1).^2 ./ var(I, 1, 1);                   % per-trial moment estimate
    fprintf('%5d   %7.3f   %7.3f   %7.2f\n', n, mean(mh), std(mh), std(mh)*sqrt(n));
end

%% 3. Figure
figure('Position', [100 100 1000 380]);
subplot(1,2,1);
Ng = logspace(0, 2, 100);
semilogx(Ng, Ng./(Ng+1), 'b-', 'LineWidth', 1.5); hold on;
semilogx(Nlist, mP, 'bo', 'MarkerFaceColor', 'b');
semilogx(Ng(Ng > 1.5), Ng(Ng > 1.5)./(Ng(Ng > 1.5)-1), 'r-', 'LineWidth', 1.5);
semilogx(Nlist(2:end), mF(2:end), 'rs', 'MarkerFaceColor', 'r');
plot([1 100], [1 1], 'k--');
ylim([0.4 2.2]); xlabel('mean scatterers per resolution cell N'); ylabel('Nakagami m');
title('Poisson count (blue) vs fixed count (red)'); grid on;
subplot(1,2,2);
semilogx(Nlist, sP, 'bo-', 'MarkerFaceColor', 'b'); hold on;
semilogx(Nlist(2:end), sF(2:end), 'rs-', 'MarkerFaceColor', 'r');
plot([1 100], sqrt(pi/(4-pi))*[1 1], 'k--');
xlabel('mean scatterers per resolution cell N'); ylabel('envelope SNR_0 = mean/std');
title('SNR_0 approaches 1.913 from both sides'); grid on; ylim([1.2 2.3]);
