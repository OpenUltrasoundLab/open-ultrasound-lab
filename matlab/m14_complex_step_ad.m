%% m14_complex_step_ad.m  --  finite differences, complex step, dual numbers, reverse mode, Wirtinger
% Purpose : 1) relative error of forward/central differences and of the complex-step
%              derivative Im f(x + i h)/h versus h, for f(x) = e^x / sqrt(sin^3 x + cos^3 x)
%              at x0 = 1.5 (a standard test function in the complex-step literature);
%           2) forward-mode AD with dual numbers [value, derivative] built from anonymous
%              functions; 3) reverse mode by hand on y = ln(x1) + x1 x2 - sin(x2) at (2, 5);
%           4) two ways the complex step silently fails (abs and the ' operator);
%           5) Wirtinger gradient dL/dw* = -A^H (y - A w) for L = ||y - A w||^2 checked
%              against finite differences, and what happens if you step along its conjugate.
% Change  : x0, the step list, the test functions.
% Expect  : best forward difference 1e-8 to 1e-9 relative error near h = 1e-8; best central
%           about 1e-11 near h = 3e-6; complex step at machine precision for every
%           h <= 1e-8 (including 1e-20 and 1e-100); dual number exact; reverse-mode
%           gradient (5.5, 1.7163); Wirtinger check about 1e-9; the conjugate step diverges.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8. Runtime about 1 s.
rng(1);

%% 1. Finite differences vs complex step
f  = @(x) exp(x) ./ sqrt(sin(x).^3 + cos(x).^3);
x0 = 1.5;
df_exact = f(x0) * (1 - (3*sin(x0)^2*cos(x0) - 3*cos(x0)^2*sin(x0)) / (2*(sin(x0)^3 + cos(x0)^3)));
fprintf('f(1.5) = %.6f, exact f''(1.5) = %.12f\n', f(x0), df_exact);
hs = 10.^(-1:-0.25:-20);
f0v = f(x0 + 0*hs);                       % same vector code path as f(x0+hs): avoids a 1-ulp mismatch
efd = abs((f(x0+hs) - f0v) ./ hs - df_exact) / abs(df_exact);
ecd = abs((f(x0+hs) - f(x0-hs)) ./ (2*hs) - df_exact) / abs(df_exact);
ecs = abs(imag(f(x0 + 1i*hs)) ./ hs - df_exact) / abs(df_exact);
pfd = max(efd, 1e-17);  pcd = max(ecd, 1e-17);  pcs = max(ecs, 1e-17);       % floor zeros for the log plot
fprintf('     h     | forward diff | central diff | complex step\n');
for h = [1e-2 1e-4 1e-6 1e-8 1e-10 1e-12 1e-16 1e-20]
    k = find(abs(log10(hs) - log10(h)) < 1e-9);
    fprintf('  %7.0e  |   %8.1e   |   %8.1e   |   %8.1e\n', h, efd(k), ecd(k), ecs(k));
end
[m1, i1] = min(efd);  [m2, i2] = min(ecd);
fprintf('best forward h = %.0e (err %.1e; sqrt(eps) = %.1e); best central h = %.0e (err %.1e; eps^(1/3) = %.1e)\n', ...
    hs(i1), m1, sqrt(eps), hs(i2), m2, eps^(1/3));
fprintf('complex step with h = 1e-100: relative error %.1e\n', abs(imag(f(x0 + 1e-100i))/1e-100 - df_exact)/df_exact);

%% 2. Forward mode with dual numbers  a = [value, derivative]
dmul = @(a, b) [a(1)*b(1), a(1)*b(2) + a(2)*b(1)];
ddiv = @(a, b) [a(1)/b(1), (a(2)*b(1) - a(1)*b(2)) / b(1)^2];
dadd = @(a, b) a + b;
dexp = @(a) [exp(a(1)), exp(a(1))*a(2)];
dsin = @(a) [sin(a(1)), cos(a(1))*a(2)];
dcos = @(a) [cos(a(1)), -sin(a(1))*a(2)];
dsqrt = @(a) [sqrt(a(1)), a(2) / (2*sqrt(a(1)))];
dcube = @(a) dmul(dmul(a, a), a);
fdual = @(a) ddiv(dexp(a), dsqrt(dadd(dcube(dsin(a)), dcube(dcos(a)))));
r = fdual([x0, 1]);                                      % seed dx/dx = 1
fprintf('\ndual number: f = %.12f, f'' = %.15f, relative error %.1e\n', r(1), r(2), abs(r(2) - df_exact)/df_exact);

%% 3. Reverse mode by hand: y = ln(x1) + x1*x2 - sin(x2) at (2, 5)
x1 = 2; x2 = 5;
v1 = log(x1); v2 = x1*x2; v3 = sin(x2); v4 = v1 + v2; y = v4 - v3;       % forward sweep (the tape)
yb = 1;                                                                   % reverse sweep: adjoints
v4b = yb;  v3b = -yb;  v1b = v4b;  v2b = v4b;
x1b = v1b/x1 + v2b*x2;  x2b = v2b*x1 + v3b*cos(x2);
fprintf('reverse mode: y = %.4f, grad = (%.4f, %.4f); analytic (1/x1 + x2, x1 - cos x2) = (%.4f, %.4f)\n', ...
    y, x1b, x2b, 1/x1 + x2, x1 - cos(x2));
fprintf('  one reverse sweep gave both partials; forward mode needs one sweep per input (2 here)\n');

%% 4. Silent failures of the complex step
g1 = @(x) x .* abs(x - 1);                % abs is not analytic: abs(z) returns |z|, a real number
fprintf('\nx|x-1| at 1.5: exact 2.0, complex step %.4f, dual %.4f\n', imag(g1(1.5 + 1e-20i))/1e-20, ...
    subsref(dmul([1.5 1], [abs(0.5), sign(0.5)*1]), struct('type', '()', 'subs', {{2}})));
g2 = @(x) x' * x;                         % ' is the CONJUGATE transpose: x'*x = |x|^2 for complex x
fprintf('x''*x at 1.5: exact 3.0, complex step %.4f;   with .'' : %.4f\n', imag(g2(1.5 + 1e-20i))/1e-20, ...
    imag((1.5 + 1e-20i).' * (1.5 + 1e-20i))/1e-20);

%% 5. Wirtinger gradient for a real loss of complex parameters
A = (randn(8,4) + 1i*randn(8,4)) / sqrt(2);  y = (randn(8,1) + 1i*randn(8,1)) / sqrt(2);
L = @(w) real((y - A*w)' * (y - A*w));
w = (randn(4,1) + 1i*randn(4,1)) / sqrt(2);
g = -A' * (y - A*w);                      % dL/dw*  (the steepest-ascent direction is 2*g in R^8)
h = 1e-6;  gfd = zeros(4,1);
for k = 1:4
    e = zeros(4,1); e(k) = h;
    gfd(k) = (L(w+e) - L(w-e))/(2*h) + 1i*(L(w+1i*e) - L(w-1i*e))/(2*h);    % dL/dRe + i dL/dIm
end
fprintf('\nWirtinger check: max |(dL/dRe + i dL/dIm) - 2 dL/dw*| / max|.| = %.1e\n', max(abs(gfd - 2*g))/max(abs(gfd)));
mu = 0.5 / max(eig(A'*A));  wa = zeros(4,1);  wb = zeros(4,1);
for it = 1:200
    wa = wa - mu * (-A' * (y - A*wa));            % correct: step along -dL/dw*
    wb = wb - mu * conj(-A' * (y - A*wb));        % bug: step along -dL/dw = -conj(dL/dw*)
end
wls = A \ y;
fprintf('after 200 steps: L(correct) = %.4f (LS optimum %.4f); L(conjugated step) = %.4g\n', L(wa), L(wls), L(wb));

%% 6. Figure
figure('Position', [100 100 640 440]);
loglog(hs, pfd, 'b.-', hs, pcd, 'r.-', hs, pcs, 'k.-'); hold on;
loglog([1e-20 1e-1], [1 1]*eps, 'k:');
set(gca, 'XDir', 'reverse'); grid on; ylim([1e-17 10]);
xlabel('step h'); ylabel('relative error of f''(1.5)');
legend('forward difference', 'central difference', 'complex step', 'machine epsilon', 'Location', 'northwest');
title('truncation error falls with h, round-off grows as eps/h; the complex step has no subtraction');
