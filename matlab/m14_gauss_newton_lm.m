%% m14_gauss_newton_lm.m  --  gradient descent vs Gauss-Newton vs Levenberg-Marquardt
% Purpose : fit a Gaussian-pulse echo model to a noisy RF echo by nonlinear least squares
%             s(t; A, t0, sg, f) = A exp(-(t-t0)^2/(2 sg^2)) cos(2 pi f (t - t0))
%           (units: t, t0, sg in us; f in MHz) and compare
%             GD + Armijo backtracking, Gauss-Newton (J'J d = -J'r), and
%             Levenberg-Marquardt ((J'J + lam diag(J'J)) d = -J'r, lam /10 on success, x10 on failure).
%           Prints the iterations, the estimates, the linearized standard deviations
%           sqrt(diag(sn^2 (J'J)^-1)) (the Gaussian-noise CRLB of Module 12 evaluated at the
%           solution), kappa(J'J) before and after column scaling, and shows that a start
%           more than half a period off in t0 converges to the wrong cycle unless the model
%           has a free carrier phase.
% Change  : theta_true, sn (noise), the starts, lam0.
% Expect  : good start: GN and LM about 10 iterations, GD + Armijo about 10^4 in raw units
%           and tens after Jacobi scaling; bad start (t0 off by 100 ns, more than T/4 = 56 ns):
%           LM lands half a period late with the polarity flipped (A < 0, t0 near 2.11 us,
%           |r| three times the noise level); the free-phase model recovers t0 = 2.00 us.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8. Runtime about 2 s.
rng(1);

%% 1. Data
fs = 50;                          % sampling rate (MHz) -> dt = 0.02 us
t = (1:1/fs:3).';                 % time axis (us), 101 samples
th_true = [1, 2.00, 0.15, 4.5];   % [A, t0 (us), sg (us), f (MHz)]
sn = 0.05;                        % noise std (peak SNR = 20 log10(1/0.05) = 26 dB)
model = @(p) p(1) * exp(-(t-p(2)).^2/(2*p(3)^2)) .* cos(2*pi*p(4)*(t-p(2)));
y = model(th_true) + sn * randn(size(t));
% analytic Jacobian of the residual r = model - y (columns: A, t0, sg, f)
jac = @(p) [exp(-(t-p(2)).^2/(2*p(3)^2)) .* cos(2*pi*p(4)*(t-p(2))), ...
     p(1)*exp(-(t-p(2)).^2/(2*p(3)^2)) .* ((t-p(2))/p(3)^2 .* cos(2*pi*p(4)*(t-p(2))) + 2*pi*p(4)*sin(2*pi*p(4)*(t-p(2)))), ...
     p(1)*exp(-(t-p(2)).^2/(2*p(3)^2)) .* cos(2*pi*p(4)*(t-p(2))) .* (t-p(2)).^2/p(3)^3, ...
    -p(1)*exp(-(t-p(2)).^2/(2*p(3)^2)) .* sin(2*pi*p(4)*(t-p(2))) .* 2*pi.*(t-p(2))];
cost = @(p) 0.5 * sum((model(p) - y).^2);
% check the Jacobian against central differences
p = [0.9 2.03 0.13 4.8];  Jn = zeros(numel(t), 4);
for k = 1:4, e = zeros(1,4); e(k) = 1e-6; Jn(:,k) = (model(p+e) - model(p-e)) / 2e-6; end
fprintf('Jacobian check: max |analytic - numeric| / max|J| = %.1e\n', max(max(abs(jac(p)-Jn))) / max(abs(Jn(:))));

%% 2. Solvers from a good start and from a start beyond a quarter period
% GD+Armijo in raw units (step capped at 0.01), GD+Armijo after Jacobi (column) scaling,
% Gauss-Newton (least-squares solve of J d = -r) and Levenberg-Marquardt.
starts = {[0.8, 2.04, 0.12, 4.8], [0.8, 2.10, 0.12, 4.5]};
labels = {'good start (dt0 = +40 ns, df = +0.3 MHz)', 'bad start (dt0 = +100 ns > T/4)'};
nm = {'GD+Armijo', 'GD+Armijo+scale', 'Gauss-Newton', 'LM'};
hist = cell(2, 4);
for s = 1:2
    fprintf('\n%s; period T = %.0f ns\n', labels{s}, 1e3/th_true(4));
    for meth = 1:4
        p = starts{s};  lam = 1e-3;  al = 1e-2;  h = cost(p);  it = 0;
        if meth == 2, Dsc = 1 ./ sqrt(sum(jac(p).^2, 1)); amax = 1; else, Dsc = ones(1,4); amax = 1e-2; end
        while it < 20000
            r = model(p) - y;  J = jac(p);  g = (J.' * r).';
            if meth <= 2                                  % (scaled) GD + Armijo backtracking
                d = -(Dsc.^2) .* g;  al = min(amax, 2*al);
                while cost(p + al*d) > cost(p) + 1e-4*al*(g*d.'), al = al/2; end
                pn = p + al*d;
            elseif meth == 3                              % Gauss-Newton
                pn = p - (J \ r).';                       % least-squares solve, J'J never formed
            else                                          % Levenberg-Marquardt
                H = J.'*J;
                while true
                    pn = p - ((H + lam*diag(diag(H))) \ g.').';
                    if cost(pn) < cost(p), lam = lam/10; break; end
                    lam = lam*10;  if lam > 1e12, pn = p; break; end
                end
            end
            it = it + 1;  step = norm(pn - p);  p = pn;  h(end+1) = cost(p);
            if step < 1e-9, break; end
        end
        hist{s, meth} = h;
        fprintf('%-16s %5d its  A = %7.4f  t0 = %.4f us  sg = %.4f us  f = %.4f MHz  |r| = %.4f\n', ...
            nm{meth}, it, p(1), p(2), abs(p(3)), abs(p(4)), sqrt(2*cost(p)));
        if s == 1 && meth == 4, p_hat = p; end
    end
end
fprintf('expected |r| for a correct fit: sn*sqrt(n) = %.3f\n', sn*sqrt(numel(t)));

%% 3. Precision at the solution and conditioning
J = jac(p_hat);  H = J.'*J;  sd = sqrt(diag(sn^2 * (H \ eye(4))));
fprintf('\nlinearized std (CRLB at the fit): A %.4f, t0 %.2f ns, sg %.2f ns, f %.4f MHz\n', sd(1), 1e3*sd(2), 1e3*sd(3), sd(4));
Dg = diag(1 ./ sqrt(diag(H)));
fprintf('kappa(J''J) = %.3g in (1, us, us, MHz) units; %.3g after column scaling\n', cond(H), cond(Dg*H*Dg));
Hs = diag([1 1e-6 1e-6 1e6]) * H * diag([1 1e-6 1e-6 1e6]);
fprintf('kappa(J''J) = %.3g if the same fit is written in SI units (1, s, s, Hz)\n', cond(Hs));

%% 4. Free carrier phase removes the cycle ambiguity in t0
modelp = @(q) q(1) * exp(-(t-q(2)).^2/(2*q(3)^2)) .* cos(2*pi*q(4)*(t-q(2)) + q(5));
costp = @(q) 0.5 * sum((modelp(q) - y).^2);
q = [starts{2} 0];  lam = 1e-3;
for it = 1:200
    r = modelp(q) - y;  Jq = zeros(numel(t), 5);
    for k = 1:5, e = zeros(1,5); e(k) = 1e-7; Jq(:,k) = (modelp(q+e) - modelp(q-e)) / 2e-7; end
    H = Jq.'*Jq;  g = Jq.'*r;
    while true
        qn = q - ((H + lam*diag(diag(H))) \ g).';
        if costp(qn) < costp(q), lam = lam/10; break; end
        lam = lam*10;  if lam > 1e12, qn = q; break; end
    end
    if norm(qn - q) < 1e-10, q = qn; break; end
    q = qn;
end
fprintf('free-phase LM from the bad start: t0 = %.4f us, phase = %.3f rad, |r| = %.4f (%d its)\n', q(2), angle(exp(1j*q(5))), sqrt(2*costp(q)), it);

%% 5. Figure
figure('Position', [100 100 1050 420]);
subplot(1,2,1);
cols = {[0.6 0.6 0.6], [0.85 0.55 0], [0 0.3 0.8], [0.8 0.2 0.5]};  hh = zeros(1,4);
for meth = 1:4
    hh(meth) = semilogy(1:numel(hist{1,meth}), hist{1,meth} - min(hist{1,4}) + 1e-12, 'Color', cols{meth}, 'LineWidth', 1.5); hold on;
end
set(gca, 'XScale', 'log'); grid on; xlim([1 2e4]);
xlabel('iteration'); ylabel('cost - final cost');
legend(hh, {'GD + Armijo (raw units)', 'GD + Armijo (Jacobi scaled)', 'Gauss-Newton', 'Levenberg-Marquardt'}, 'Location', 'southwest');
title('good start: linear vs quadratic convergence');
subplot(1,2,2);
plot(t, y, 'Color', [0.6 0.6 0.6]); hold on;
plot(t, model(p_hat), 'b', 'LineWidth', 1.5);
plot(t, model(starts{2}), 'r--');
xlabel('time (\mus)'); ylabel('amplitude'); grid on; xlim([1.4 2.6]);
legend('data', 'LM fit (good start)', 'bad start model', 'Location', 'southwest');
title('echo, fit, and a start that lands half a period late');
