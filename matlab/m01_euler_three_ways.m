%% m01_euler_three_ways.m  --  e^{j theta} = cos(theta) + j sin(theta), three ways
% Purpose : verify Euler's formula numerically by the three routes of section 01.3:
%   (1) Taylor series: partial sums of sum (j theta)^k / k!
%   (2) ODE: f' = j f means velocity is position turned by 90 deg, so |f| stays 1
%   (3) limit: (1 + j theta/n)^n = n small rotate-and-scale steps
%           (this is exactly forward-Euler integration of f' = j f with step theta/n)
% Change  : theta (try pi, 2*pi, 10); the list of n.
% Expect  : theta = 1: 8 terms give error 2.5e-5; theta = pi: 13 terms for error < 1e-3.
%           Limit at theta = pi: n = 4 -> |.| = 2.614, angle 152.6 deg;
%           n = 100 -> |.| = 1.0506 (close to exp(theta^2/(2n)) = 1.0506), angle 179.94 deg.
% Runtime : < 1 s. Base MATLAB / GNU Octave, no toolboxes.
clear; close all; rng(1);

%% 1. Taylor series, term by term
for theta = [1, pi]
    exact = exp(1j*theta);
    S = 0;
    fprintf('Taylor, theta = %.4f rad   (exact %.4f %+.4fj)\n', theta, real(exact), imag(exact));
    fprintf('  terms   term (j theta)^k/k!      partial sum          |error|\n');
    for k = 0:13
        term = (1j*theta)^k / factorial(k);   % powers of j cycle 1, j, -1, -j
        S = S + term;
        fprintf('  %4d   %+9.4f %+9.4fj   %+8.4f %+8.4fj   %9.2e\n', k+1, real(term), imag(term), real(S), imag(S), abs(S - exact));
    end
end
fprintf('\n');

%% 2. ODE: velocity is perpendicular to position
theta = linspace(0, 2*pi, 2001);
f  = exp(1j*theta);
h  = theta(2) - theta(1);
df = (f(3:end) - f(1:end-2)) / (2*h);         % central difference d f / d theta
fm = f(2:end-1);
fprintf('ODE check: max |df/dtheta - j f| = %.2e (finite-difference error ~ h^2/6)\n', max(abs(df - 1j*fm)));
fprintf('           max |Re(conj(f) df)| = %.2e  (velocity _|_ position, so |f| is constant)\n', max(abs(real(conj(fm).*df))));
fprintf('           max ||f| - 1| = %.2e ; speed |df| = %.6f rad per rad\n\n', max(abs(abs(f) - 1)), mean(abs(df)));

%% 3. Limit (1 + j theta / n)^n and forward Euler
theta = pi;
ns = [1 2 4 8 16 32 100 1000 10000];
lim_err = zeros(size(ns));
fprintf('Limit, theta = pi\n     n     |(1+j pi/n)^n|   exp(pi^2/2n)   angle (deg)   |error|\n');
for i = 1:numel(ns)
    n = ns(i);
    v = (1 + 1j*theta/n)^n;
    y = 1; for k = 1:n, y = y + (theta/n)*1j*y; end   % forward Euler on f' = j f
    lim_err(i) = abs(v - exp(1j*theta));
    fprintf('%6d   %12.4f   %12.4f   %10.2f   %9.2e   (Euler-step diff %.1e)\n', n, abs(v), exp(theta^2/(2*n)), angle(v)*180/pi, lim_err(i), abs(y - v));
end
fprintf('(1 + j/2)^2 by hand = 1 + j - 1/4 = %.4f %+.4fj\n', real((1 + 0.5j)^2), imag((1 + 0.5j)^2));

%% 4. Figures
uc = exp(1j*linspace(0, 2*pi, 400));
figure('Position', [80 80 1100 380]);
subplot(1,3,1); hold on; grid on;
plot(real(uc), imag(uc), ':', 'Color', [0.5 0.5 0.5]);
S = cumsum((1j*pi).^(0:14) ./ factorial(0:14));
plot(real([0 S]), imag([0 S]), '-o', 'Color', [0.12 0.37 0.82], 'MarkerSize', 4, 'LineWidth', 1.5);
plot(-1, 0, 'p', 'MarkerSize', 12, 'Color', [0.76 0.21 0.44]);
axis equal; axis([-4.5 1.5 -2.5 3.5]); xlabel('real part'); ylabel('imaginary part');
title('Taylor partial sums, \theta = \pi  ->  -1');
subplot(1,3,2); hold on; grid on;
plot(real(uc), imag(uc), ':', 'Color', [0.5 0.5 0.5]);
cols = [0.76 0.21 0.44; 0.66 0.43 0.02; 0.12 0.37 0.82];
nlist = [4 8 32];
for i = 1:3
    n = nlist(i); pts = (1 + 1j*pi/n).^(0:n);
    plot(real(pts), imag(pts), '-o', 'Color', cols(i,:), 'MarkerSize', 3, 'LineWidth', 1.3);
end
axis equal; axis([-2.6 1.4 -0.6 3]); xlabel('real part'); ylabel('imaginary part');
legend('unit circle', 'n = 4', 'n = 8', 'n = 32', 'Location', 'northwest');
title('(1 + j\pi/n)^k, k = 0..n');
subplot(1,3,3);
Nt = 1:20; terr = zeros(size(Nt));
for i = Nt, terr(i) = abs(sum((1j*pi).^(0:i-1) ./ factorial(0:i-1)) - (-1)); end
loglog(Nt, terr, 'o-', 'Color', [0.12 0.37 0.82], 'LineWidth', 1.5); hold on; grid on;
loglog(ns, lim_err, 's-', 'Color', [0.76 0.21 0.44], 'LineWidth', 1.5);
loglog(ns, pi^2./(2*ns), 'k--');
xlabel('number of terms N  or  steps n'); ylabel('|error| at \theta = \pi');
legend('Taylor, N terms', 'limit, n steps', '\pi^2/(2n)', 'Location', 'southwest');
title('Taylor converges factorially, the limit like 1/n');
