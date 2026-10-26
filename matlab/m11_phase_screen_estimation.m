%% m11_phase_screen_estimation.m
% Purpose : near-field phase-screen (time-delay) estimation from diffuse
%           speckle on a 64-element aperture. Each element sees partially
%           correlated speckle (adjacent correlation rho) delayed by a smooth
%           random screen tau_i. Adjacent-element NCC delays are integrated
%           with a cumulative sum (Flax & O'Donnell style) and, for
%           comparison, with a least-squares fit to lags 1..3.
% Change  : rho, snrdB, screen rms and correlation length, nLagLS.
% Expect  : rho = 0.984, 20 dB: adjacent-pair jitter sigma_D = 1.7 ns,
%           close to the decorrelation-only Walker-Trahey floor (1.6 ns for
%           1/rho^2 - 1 = 0.033, T = 2 us). The cumulative-sum error grows as
%           sqrt(n)*sigma_D (a random walk): 12 ns rms (22 deg at 5 MHz) at
%           element 64 against sqrt(63)*sigma_D = 13.8 ns. Weighted least
%           squares over lags 1..3 improves this by only about 5 %: with the
%           Markov (rho^|lag|) inter-element correlation used here, a lag-l
%           pair is as noisy as l adjacent pairs and adds little independent
%           information. With rho = 0.95, about 1 % of lag-3 pairs hop by
%           1/f0 and the LS result (about 49 ns) is worse than the plain
%           cumulative sum (23 ns). When electronic noise dominates instead
%           (rho = 0.9999, snrdB = 10) LS over lags 1..3 wins clearly
%           (about 3 ns vs 7 ns). Noise is white over the full 0-20 MHz band
%           (no band-pass prefilter), which inflates sigma_D at low SNR.
rng(1);

%% Parameters
f0 = 5e6;  fs = 40e6;  B6 = 0.6;   % centre frequency, sampling (Hz), -6 dB fractional BW
Nel = 64;                          % elements (pitch 0.3 mm)
N  = 1024;  n0 = 400;  W = 80;     % record, kernel start (1-based), kernel 2 us
K  = 24;                           % search +-K lags (+-600 ns)
rho = 0.984;                       % adjacent-element speckle correlation (VCZ at focus, uniform 64-el. transmit: 1 - 1/64)
snrdB = 20;                        % channel SNR (dB, real-part power ratio over the record)
rmsNs = 30;  Lcorr = 6;            % screen rms (ns) and correlation length (elements)
nReal = 100;                       % realizations for the error statistics
nLagLS = 3;                        % lags used by the least-squares integrator

f  = (0:N-1)' * fs / N;  f(f >= fs/2) = f(f >= fs/2) - fs;          % signed frequencies (Hz)
sf = B6 * f0 / 2 / sqrt(2*log(2));
H  = exp(-(abs(f) - f0).^2 / (2*sf^2));  H(N/2+1) = 0;
fitCos = @(Rm, R0, Rp) atan((Rp - Rm) ./ (2*R0 .* sin(acos((Rm + Rp) ./ (2*R0))))) ./ acos((Rm + Rp) ./ (2*R0));
idxW = (n0:n0+W-1)';  lags = -K:K;  Yidx = bsxfun(@plus, idxW, lags);
el = (0:Nel-1)';
% Gaussian-correlated screen: filter white noise along the aperture
gk = exp(-(-3*Lcorr:3*Lcorr).^2 / (2*(Lcorr/2)^2));  gk = gk / sqrt(sum(gk.^2));
% Least-squares design matrix: tau(i+l) - tau(i) = D(i,l), tau(1) = 0
rows = [];  for l = 1:nLagLS, for i = 1:Nel-l, rows = [rows; i, i+l, l]; end, end %#ok<AGROW>
A = zeros(size(rows,1), Nel);
for r = 1:size(rows,1), A(r, rows(r,1)) = -1; A(r, rows(r,2)) = 1; end
A = A(:, 2:end);                                                     % fix tau(1) = 0

errCS = zeros(Nel, nReal);  errLS = errCS;  dAdj = zeros(Nel-1, nReal);  trueAdj = dAdj;
for r = 1:nReal
    tau = conv(randn(Nel + 6*Lcorr, 1), gk', 'valid');  tau = tau(1:Nel);
    tau = (tau - mean(tau)) / std(tau) * rmsNs * 1e-9;               % screen (s)
    % partially correlated speckle across elements: AR(1) in element index
    S = zeros(N, Nel);  S(:,1) = randn(N,1);
    for i = 2:Nel, S(:,i) = rho*S(:,i-1) + sqrt(1 - rho^2)*randn(N,1); end
    X = real(ifft(bsxfun(@times, fft(S), H) .* exp(-1j*2*pi*f*tau')));  % delayed band-limited channels
    X = X / sqrt(mean(X(:).^2));
    X = X + randn(N, Nel) * 10^(-snrdB/20);
    D = zeros(Nel, nLagLS);  Rpk = zeros(Nel, nLagLS);
    for l = 1:nLagLS
        for i = 1:Nel-l
            xw = X(idxW, i);  Yw = X(:, i+l);  Yw = Yw(Yidx);
            R = (xw' * Yw) ./ sqrt(sum(xw.^2) * sum(Yw.^2, 1));
            [~, im] = max(R(2:end-1));  im = im + 1;
            dl = fitCos(R(im-1), R(im), R(im+1));  if ~isfinite(dl), dl = 0; end
            D(i, l) = (lags(im) + dl) / fs;                         % estimated tau(i+l) - tau(i)
            Rpk(i, l) = R(im);                                      % peak NCC (quality)
        end
    end
    tauCS = [0; cumsum(D(1:Nel-1, 1))];                              % cumulative sum of adjacent delays
    % weighted least squares: weight = Fisher-information proxy rho^2/(1-rho^2) from the peak NCC,
    % then one pass of outlier rejection (residual > 1/(4 f0)) to drop side-peak hops
    b = zeros(size(rows,1), 1);  wq = b;
    for q = 1:size(rows,1), b(q) = D(rows(q,1), rows(q,3)); rq = min(Rpk(rows(q,1), rows(q,3)), 0.999); wq(q) = rq^2/(1 - rq^2); end
    sw = sqrt(wq);  t1 = (bsxfun(@times, A, sw)) \ (b .* sw);
    keep = abs(A*t1 - b) < 0.25/f0;
    if rank(A(keep,:)) < Nel-1, keep = true(size(b)); end          % never drop an element entirely
    t1 = (bsxfun(@times, A(keep,:), sw(keep))) \ (b(keep) .* sw(keep));
    tauLS = [0; t1];
    errCS(:, r) = tauCS - (tau - tau(1));
    errLS(:, r) = tauLS - (tau - tau(1));
    dAdj(:, r) = D(1:Nel-1, 1);  trueAdj(:, r) = diff(tau);
end

%% Statistics
eAdj = dAdj(:) - trueAdj(:);
hop = abs(eAdj) > 0.5/f0;
sigD = std(eAdj(~hop));
sdCS = sqrt(mean(errCS.^2, 2));  sdLS = sqrt(mean(errLS.^2, 2));
fprintf('rho = %.3f, SNR = %d dB, screen rms = %d ns, Lcorr = %.0f el, %d realizations\n', rho, snrdB, rmsNs, Lcorr, nReal);
fprintf('adjacent-pair jitter sigma_D = %.2f ns (peak hops: %.2f %% of pairs)\n', sigD*1e9, 100*mean(hop));
fprintf('rms error at element 64: cumsum = %.1f ns (random walk sqrt(63)*sigma_D = %.1f ns), weighted LS lags 1..%d = %.1f ns\n', ...
    sdCS(end)*1e9, sqrt(63)*sigD*1e9, nLagLS, sdLS(end)*1e9);
fprintf('phase error at element 64, cumsum: %.0f deg at %.0f MHz\n', sdCS(end)*f0*360, f0/1e6);

%% Figure
figure;
subplot(1, 2, 1);
plot(el, (tau - tau(1))*1e9, 'k-', el, tauCS*1e9, 'o-', el, tauLS*1e9, 's-');
xlabel('element'); ylabel('delay re element 0 (ns)'); grid on;
legend('true screen', 'cumsum', sprintf('WLS, lags 1..%d', nLagLS), 'Location', 'southwest');
title('Last realization');
subplot(1, 2, 2);
plot(el, sdCS*1e9, 'o-', el, sdLS*1e9, 's-', el, sqrt(el)*sigD*1e9, 'k--');
xlabel('element n'); ylabel('rms error (ns)'); grid on;
legend('cumsum', 'weighted LS', 'sqrt(n) \sigma_D', 'Location', 'northwest');
title('Error growth');
