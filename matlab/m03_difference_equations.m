%% m03_difference_equations.m
% Purpose : run the three basic recursions (moving average, first
%           difference, one-pole smoother) by explicit loops and with
%           filter, print the by-hand samples, and check the standard
%           numbers: -3 dB points, group delay, noise variance gain,
%           moving-average sidelobes, instability, pole-pair geometry.
% Change  : a (one-pole coefficient), L (moving-average length).
% Expect  : step response 0.2 0.36 0.488 0.5904 0.67232; a = 0.9: -3 dB at
%           0.1055 rad/sample, group delay 9.00 samples at DC, noise gain
%           0.0526; a = 1.05 doubles every 14.21 samples; the pole pair
%           0.9 exp(+-j pi/4) peaks 12.04 dB above DC.
clear; close all; rng(1);

%% One-pole smoother y(n) = a y(n-1) + (1-a) x(n), by loop and by filter
a = 0.8; N = 40;
x = ones(1, N);                              % step input
y = zeros(1, N); yprev = 0;
for n = 1:N
    y(n) = a * yprev + (1 - a) * x(n);
    yprev = y(n);
end
yf = filter(1 - a, [1 -a], x);
fprintf('one-pole a = 0.8, step response y(0..5): %s\n', mat2str(round(y(1:6)*1e5)/1e5));
fprintf('  closed form 1 - a^(n+1)            : %s\n', mat2str(round((1 - a.^(1:6))*1e5)/1e5));
fprintf('  loop vs filter max diff: %.1e; time constant -1/ln(a) = %.2f samples\n', max(abs(y - yf)), -1/log(a));
hi = filter(0.5, [1 -0.5], [1 zeros(1, 5)]);
fprintf('one-pole a = 0.5, impulse response: %s\n', mat2str(hi));

%% Moving average and first difference, by loop and by filter
xr = [zeros(1,3) ones(1,4) zeros(1,4)];      % a rect
ma = zeros(size(xr)); fd = zeros(size(xr));
xa = @(k) (k >= 1) * xr(max(k, 1));          % x with zeros before the start
for n = 1:numel(xr)
    ma(n) = (xa(n) + xa(n-1) + xa(n-2)) / 3;
    fd(n) = xa(n) - xa(n-1);
end
fprintf('rect                 : %s\n', mat2str(xr));
fprintf('3-pt average (loop)  : %s\n', mat2str(round(ma*1e4)/1e4));
fprintf('3-pt average (filter): %s\n', mat2str(round(filter(ones(1,3)/3, 1, xr)*1e4)/1e4));
fprintf('first difference     : %s  (edges only)\n', mat2str(fd));
fprintf('conv([1 -1], [1 1 1 1]) = %s\n', mat2str(conv([1 -1], [1 1 1 1])));

%% Standard numbers for the one-pole smoother
wg = linspace(0, pi, 200001);
for a = [0.8 0.9]
    H = (1 - a) ./ (1 - a * exp(-1j * wg));
    w3 = acos(1 - (1 - a)^2 / (2 * a));          % closed form
    [~, i3] = min(abs(abs(H) - 1/sqrt(2)));      % numerical check
    ph = unwrap(angle(H)); gd0 = -(ph(2) - ph(1)) / (wg(2) - wg(1));
    hh = (1 - a) * a.^(0:2000);
    fprintf('a = %.1f: -3 dB at %.4f (grid %.4f) rad/sample; |H(pi)| = %.4f (%.1f dB); gd(0) = %.2f = a/(1-a) = %.2f\n', ...
            a, w3, wg(i3), abs(H(end)), 20*log10(abs(H(end))), gd0, a/(1-a));
    fprintf('        noise variance gain sum(h^2) = %.4f = (1-a)/(1+a) = %.4f; tau = %.2f samples\n', ...
            sum(hh.^2), (1-a)/(1+a), -1/log(a));
end

%% Moving-average sidelobes and the first-difference approximation
for L = [3 5 8]
    Hm = abs(sin(L*wg/2) ./ (L*sin(wg/2))); Hm(1) = 1;
    fprintf('L = %d moving average: first sidelobe %.2f dB, -3 dB at %.3f rad/sample\n', ...
            L, 20*log10(max(Hm(wg > 2*pi/L))), wg(find(Hm < 1/sqrt(2), 1)));
end
fprintf('first difference at w = pi/2: |H| = %.4f, angle = %.1f deg, 2 sin(w/2)/w = %.3f\n', ...
        abs(1 - exp(-1j*pi/2)), angle(1 - exp(-1j*pi/2))*180/pi, 2*sin(pi/4)/(pi/2));

%% Instability: pole at 1.05
hu = filter(1, [1 -1.05], [1 zeros(1, 40)]);
fprintf('a = 1.05: h(40) = %.2f, doubling time ln2/ln1.05 = %.2f samples\n', hu(41), log(2)/log(1.05));

%% Pole pair 0.9 exp(+-j pi/4): |H| from distances
p = 0.9 * exp(1j*pi/4);
A = real(poly([p conj(p)]));                     % [1, -2 r cos(theta), r^2]
Hp = 1 ./ polyval(A, exp(1j*wg)) .* exp(1j*wg).^2;   % H(z) = 1/(1 + a1 z^-1 + a2 z^-2)
e = exp(1j*pi/4);
fprintf('pole pair: distances at w = pi/4: %.3f, %.3f -> |H| = %.3f; at DC: %.3f -> |H| = %.3f\n', ...
        abs(e - p), abs(e - conj(p)), 1/(abs(e-p)*abs(e-conj(p))), abs(1-p), 1/abs(1-p)^2);
[pk, ip] = max(abs(Hp));
fprintf('  peak %.3f at w = %.3f pi, %.2f dB above DC; sum|h| = %.2f\n', pk, wg(ip)/pi, 20*log10(pk/abs(Hp(1))), ...
        sum(abs(filter(1, A, [1 zeros(1, 2000)]))));

%% Figures
figure('Name', 'Difference equations');
subplot(2,2,1); stem(0:39, yf, 'filled'); xlabel('n (samples)'); title('one-pole a = 0.8: step response');
subplot(2,2,2); stem(0:40, hu, 'r', 'filled'); hold on; stem(0:40, 0.9.^(0:40), 'b'); hold off;
xlabel('n (samples)'); title('p^n: 0.9 (b) decays, 1.05 (r) grows');
subplot(2,2,3);
H8 = 0.2 ./ (1 - 0.8*exp(-1j*wg)); H9 = 0.1 ./ (1 - 0.9*exp(-1j*wg));
plot(wg, abs(H8), 'b', wg, abs(H9), 'r', [0 pi], [1 1]/sqrt(2), 'k:'); xlim([0 1]);
xlabel('\omega (rad/sample)'); title('|H|: a = 0.8 (b), 0.9 (r)');
subplot(2,2,4); plot(wg/pi, 20*log10(abs(Hp)/abs(Hp(1)))); grid on;
xlabel('\omega/\pi'); ylabel('dB re DC'); title('pole pair 0.9e^{\pm j\pi/4}');
