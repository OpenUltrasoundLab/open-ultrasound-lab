%% m36_review_clock.m  510(k) calendar time = FDA days + holds (Monte Carlo, illustrative)
% Purpose : Equation (36.5): T_cal = T_FDA + sum(T_hold). FDA's MDUFA V goals concern FDA days
%           (decision within 90 FDA days for 95 % of 510(k)s); holds stop the clock. A hold must be
%           resolved within 180 calendar days or the 510(k) is considered withdrawn.
%           This script draws many hypothetical submissions under STATED, ILLUSTRATIVE assumptions
%           and shows how the sponsor's response time drives calendar time and withdrawals.
% Assumptions (illustrative inputs, NOT FDA statistics; change them):
%   pScreen  probability the eSTAR fails technical screening; fix time uniform 5-30 days
%   pAI      probability of an AI (additional information) letter; issued on FDA day 40-60
%   response time lognormal with median medResp(m) and log-sd sigResp
%   FDA days after the response: uniform 15-30; without an AI letter: decision on FDA day 45-90
% Expected: median and 90th percentile of calendar days (withdrawn = never decided; p90 prints
%           'withdrawn' when more than 10 % are withdrawn) rise with the response median;
%           withdrawals appear when the response tail passes 180 days.
% Base MATLAB / GNU Octave only. Runtime ~1 s.
clear; close all; rng(2);

N = 20000;
pScreen = 0.2; pAI = 0.6; sigResp = 0.5;
medResp = [45 90 150];                      % three response models, days
names = {'fast response', 'typical response', 'slow response'};
pct = @(v, q) v(max(1, min(numel(v), round(q*numel(v)))));   % percentile of a sorted vector

figure('Color', 'w');
fprintf('Illustrative 510(k) clock model, N = %d per case, pAI = %.2f, pScreen = %.2f\n', N, pAI, pScreen);
fprintf('%-18s %10s %10s %10s %10s %10s\n', 'case', 'FDA med', 'cal med', 'cal p90', 'withdrawn', '>90 FDA d');
for m = 1:numel(medResp)
  scr  = rand(N, 1) < pScreen;  tScr = scr .* (5 + 25*rand(N, 1));
  ai   = rand(N, 1) < pAI;
  tAIday = 40 + 20*rand(N, 1);
  tResp  = exp(log(medResp(m)) + sigResp*randn(N, 1));
  tAfter = 15 + 15*rand(N, 1);
  tNoAI  = 45 + 45*rand(N, 1);
  tFDA = ai .* (tAIday + tAfter) + (~ai) .* tNoAI;
  tHold = tScr + ai .* tResp;
  wd = (ai & tResp > 180);                 % response later than 180 days -> withdrawn
  tCal = tFDA + tHold;
  tAll = tCal; tAll(wd) = Inf;             % withdrawn: never decided (counted, not dropped)
  ok = sort(tAll); f = sort(tFDA(~wd));
  p90 = pct(ok, 0.9);
  if isinf(p90), s90 = 'withdrawn'; else, s90 = sprintf('%.0f', p90); end   % > 10 % withdrawn
  fprintf('%-18s %10.0f %10.0f %10s %9.1f%% %9.1f%%\n', names{m}, pct(f, 0.5), pct(ok, 0.5), s90, 100*mean(wd), 100*mean(tFDA(~wd) > 90));
  subplot(numel(medResp), 1, m);
  hist(tCal(~wd), 0:10:400); xlim([0 400]);
  ylabel('count'); title(sprintf('%s (median response %d d)', names{m}, medResp(m)));
  if m == numel(medResp), xlabel('calendar days from receipt to decision'); end
end
fprintf('\nThe FDA-day distribution is the same in every case; only the hold changes.\n');
