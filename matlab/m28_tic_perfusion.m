%% m28_tic_perfusion.m  --  Time-intensity curves: bolus (lognormal) and flash-replenishment fits
% Purpose : simulate ROI-averaged linear intensity with speckle (gamma-distributed, shape Ns)
%           and fit (1) a lognormal bolus model and (2) the Wei flash-replenishment model
%           with a hand-written Levenberg-Marquardt loop (numeric Jacobian). Then show why
%           fitting log-compressed display values biases the perfusion parameters.
% Models  : bolus  I(t) = AUC / ((t - t0) s sqrt(2 pi)) exp(-(ln(t - t0) - m)^2 / (2 s^2)) + bg
%           MTT = exp(m + s^2/2) (after t0), TTP = t0 + exp(m - s^2), PE = peak above bg
%           flash  I(t) = A (1 - exp(-beta t)) + bg, flow index A * beta
% Change  : Ns (speckle cells in the ROI), frame rate, true parameters, DR of the display map.
% Expect  : linear-domain fits recover MTT = 11.97 s and TTP = 9.98 s within a few percent
%           (Monte Carlo std printed); fitting the 40 dB log-compressed values lengthens MTT
%           by a factor of about 2 and the flash fit recovers A, beta to a few percent.
% Provenance: Monte Carlo speckle on model curves (lognormal, Wei), not patient data.
% Runtime : about 2 s.
clear; close all; rng(1);

%% Bolus truth and acquisition
fr = 4;  t = (0:1/fr:60)';                        % s, 4 frames per second, 60 s
th = [1000, log(10), 0.6, 3, 2];                  % [AUC (a.u. s), m, s, t0 (s), bg]
lgn = @(p, t) (t > p(4)) .* p(1) ./ (max(t - p(4), 1e-9) * p(3) * sqrt(2*pi)) ...
              .* exp(-(log(max(t - p(4), 1e-9)) - p(2)).^2 / (2 * p(3)^2)) + p(5);
Ns = 50;                                          % independent speckle cells averaged in the ROI
noisy = @(I) I .* (-sum(log(rand(numel(I), Ns)), 2) / Ns);   % mean I, std I/sqrt(Ns)
mtt = @(p) exp(p(2) + p(3)^2/2); ttp = @(p) p(4) + exp(p(2) - p(3)^2);
fprintf('Truth: MTT = %.2f s, TTP = %.2f s, peak above bg = %.1f, AUC = %.0f\n', mtt(th), ttp(th), max(lgn(th, t)) - th(5), th(1));

%% Levenberg-Marquardt, written inline (scripts cannot define local functions portably)
keepFit = cell(1, 2);
nTrial = 100; est = zeros(nTrial, 4); estLog = zeros(nTrial, 2);
DR = 40;                                          % display dynamic range (dB)
for trial = 1:nTrial
  y = noisy(lgn(th, t));
  for mode = 1:2
    if mode == 1
      yy = y;
    else                                          % naive: fit the log-compressed gray levels
      Iref = max(y); yy = 255 * min(max((10*log10(y / Iref) + DR) / DR, 0), 1);
    end
    [pk, ip] = max(yy); base = median(yy(1:10));
    p = [ (pk - base) * 8, log(max(t(ip) - 3, 1)), 0.6, 3, base ];   % initial guess
    lam = 1e-2; r = yy - lgn(p, t); cost = r' * r;
    for it = 1:200
      J = zeros(numel(t), 5);
      for k = 1:5
        dp = zeros(1, 5); dp(k) = 1e-6 * max(1, abs(p(k)));
        J(:, k) = (lgn(p + dp, t) - lgn(p - dp, t)) / (2 * dp(k));
      end
      Hm = J' * J; g = J' * r;
      stepOK = false;
      while ~stepOK && lam < 1e12
        d = (Hm + lam * diag(diag(Hm) + 1e-12)) \ g;
        pn = p + d'; pn(3) = max(pn(3), 0.05); pn(4) = min(pn(4), t(ip) - 0.5);
        rn = yy - lgn(pn, t); cn = rn' * rn;
        if cn < cost, p = pn; r = rn; dc = cost - cn; cost = cn; lam = lam / 3; stepOK = true; else, lam = lam * 4; end
      end
      if ~stepOK || dc < 1e-10 * cost, break; end
    end
    if mode == 1, est(trial, :) = [mtt(p), ttp(p), max(lgn(p, t)) - p(5), p(1)]; else, estLog(trial, :) = [mtt(p), ttp(p)]; end
    if trial == 1, keepFit{mode} = struct('y', yy, 'p', p); end
  end
end
fprintf('Linear-intensity fits (%d trials, Ns = %d): MTT %.2f +- %.2f s, TTP %.2f +- %.2f s, peak %.1f +- %.1f, AUC %.0f +- %.0f\n', ...
  nTrial, Ns, mean(est(:,1)), std(est(:,1)), mean(est(:,2)), std(est(:,2)), mean(est(:,3)), std(est(:,3)), mean(est(:,4)), std(est(:,4)));
fprintf('Log-compressed fits (DR = %d dB):        MTT %.2f +- %.2f s (x%.2f), TTP %.2f +- %.2f s\n', ...
  DR, mean(estLog(:,1)), std(estLog(:,1)), mean(estLog(:,1)) / mtt(th), mean(estLog(:,2)), std(estLog(:,2)));

%% Flash-replenishment (Wei et al. 1998 model)
tf = (0:1/fr:15)'; thF = [50, 0.8, 2];            % A, beta (1/s), bg
mdlF = @(p, t) p(1) * (1 - exp(-p(2) * t)) + p(3);
estF = zeros(nTrial, 3);
for trial = 1:nTrial
  y = noisy(mdlF(thF, tf));
  p = [max(y) - y(1), 1, y(1)]; lam = 1e-2; r = y - mdlF(p, tf); cost = r' * r;
  for it = 1:200
    J = zeros(numel(tf), 3);
    for k = 1:3
      dp = zeros(1, 3); dp(k) = 1e-6 * max(1, abs(p(k)));
      J(:, k) = (mdlF(p + dp, tf) - mdlF(p - dp, tf)) / (2 * dp(k));
    end
    Hm = J' * J; g = J' * r; stepOK = false;
    while ~stepOK && lam < 1e12
      d = (Hm + lam * diag(diag(Hm) + 1e-12)) \ g; pn = p + d'; pn(2) = max(pn(2), 1e-3);
      rn = y - mdlF(pn, tf); cn = rn' * rn;
      if cn < cost, p = pn; r = rn; dc = cost - cn; cost = cn; lam = lam / 3; stepOK = true; else, lam = lam * 4; end
    end
    if ~stepOK || dc < 1e-10 * cost, break; end
  end
  estF(trial, :) = p;
  if trial == 1, keepF = struct('y', y, 'p', p); end
end
fprintf('Flash-replenishment fits: A %.1f +- %.1f (true 50), beta %.3f +- %.3f 1/s (true 0.8), A*beta %.1f +- %.1f (true 40)\n', ...
  mean(estF(:,1)), std(estF(:,1)), mean(estF(:,2)), std(estF(:,2)), mean(estF(:,1).*estF(:,2)), std(estF(:,1).*estF(:,2)));

%% Figures
figure('Position', [50 50 960 330]);
subplot(1, 3, 1); plot(t, keepFit{1}.y, '.', 'Color', [0.6 0.6 0.6]); hold on; plot(t, lgn(th, t), 'k--', t, lgn(keepFit{1}.p, t), 'b'); grid on;
xlabel('time (s)'); ylabel('linear intensity (a.u.)'); title('bolus: lognormal fit'); legend('data', 'truth', 'fit');
subplot(1, 3, 2); plot(t, keepFit{2}.y, '.', 'Color', [0.6 0.6 0.6]); hold on; plot(t, lgn(keepFit{2}.p, t), 'r'); grid on;
xlabel('time (s)'); ylabel('gray level (40 dB map)'); title('same data, log-compressed');
subplot(1, 3, 3); plot(tf, keepF.y, '.', 'Color', [0.6 0.6 0.6]); hold on; plot(tf, mdlF(thF, tf), 'k--', tf, mdlF(keepF.p, tf), 'b'); grid on;
xlabel('time after flash (s)'); ylabel('linear intensity (a.u.)'); title('flash-replenishment');
