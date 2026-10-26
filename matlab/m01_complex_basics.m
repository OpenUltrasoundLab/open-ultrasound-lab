%% m01_complex_basics.m  --  complex numbers as arrows
% Purpose : check by computer every hand calculation of Module 01, sections 01.1-01.2:
%           rectangular <-> polar, addition (tip to tail), multiplication
%           (rotate and scale), conjugate (mirror), |z|^2 = z z*, division,
%           powers of j, and the atan vs atan2 pitfall.
% Change  : z and w below; try w = 1j (pure quarter turn), w = -1 (half turn),
%           w = 2 (pure scale), w = exp(1j*pi/6) (pure 30 deg rotation).
% Expect  : |z| = 5, angle 53.13 deg; z*w = -5 + 10j, |zw| = 11.180 = 5 x 2.236,
%           angle 116.57 = 53.13 + 63.43 deg; z/w = 2.2 - 0.4j; z z* = 25;
%           atan(4/3) = 53.13 deg but the true angle of -3-4j is -126.87 deg.
% Runtime : < 1 s. Base MATLAB / GNU Octave, no toolboxes.
clear; close all; rng(1);

%% 1. One number, two descriptions
z = 3 + 4j;                       % rectangular: 3 units right, 4 units up
w = 1 + 2j;
r_z  = abs(z);                    % length  sqrt(3^2 + 4^2)
th_z = angle(z);                  % angle   atan2(4, 3), in rad, range (-pi, pi]
fprintf('z = %g %+gj  ->  |z| = %.4f, angle = %.4f rad = %.2f deg\n', real(z), imag(z), r_z, th_z, th_z*180/pi);
fprintf('w = %g %+gj  ->  |w| = %.4f, angle = %.2f deg\n', real(w), imag(w), abs(w), angle(w)*180/pi);
z_back = r_z*cos(th_z) + 1j*r_z*sin(th_z);   % polar -> rectangular
fprintf('back to rectangular: %.4f %+.4fj (Euler: r e^{j theta})\n\n', real(z_back), imag(z_back));

%% 2. Arithmetic by hand vs built-in
a = real(z); b = imag(z); c = real(w); d = imag(w);
s_hand = (a + c) + 1j*(b + d);                 % add components (tip to tail)
p_hand = (a*c - b*d) + 1j*(a*d + b*c);         % FOIL with j^2 = -1
q_hand = ((a*c + b*d) + 1j*(b*c - a*d)) / (c^2 + d^2);   % z w* / |w|^2
fprintf('z + w = %g %+gj   (built-in %g %+gj)\n', real(s_hand), imag(s_hand), real(z+w), imag(z+w));
fprintf('z * w = %g %+gj   (built-in %g %+gj)\n', real(p_hand), imag(p_hand), real(z*w), imag(z*w));
fprintf('   |zw| = %.4f = |z||w| = %.4f ; angle(zw) = %.2f = %.2f + %.2f deg\n', ...
        abs(z*w), abs(z)*abs(w), angle(z*w)*180/pi, angle(z)*180/pi, angle(w)*180/pi);
fprintf('z / w = %g %+gj   (built-in %g %+gj)\n', real(q_hand), imag(q_hand), real(z/w), imag(z/w));
fprintf('   |z/w| = %.4f = |z|/|w| ; angle = %.2f deg = 53.13 - 63.43\n', abs(z/w), angle(z/w)*180/pi);
fprintf('conj(z) = %g %+gj (mirror in the real axis);  z*conj(z) = %g %+gj = |z|^2\n', ...
        real(conj(z)), imag(conj(z)), real(z*conj(z)), imag(z*conj(z)));
fprintf('but z^2 = %g %+gj  (squaring is not |z|^2)\n\n', real(z^2), imag(z^2));

%% 3. j is the quarter-turn operator
for k = 0:4
    v = [round(real(1j^k)), round(imag(1j^k))];
    v(v == 0) = 0;                % print +0, not -0
    fprintf('j^%d = %+g %+gj  (angle %4.0f deg)\n', k, v(1), v(2), k*90);
end
fprintf('1/j = %g %+gj  (a quarter turn backwards)\n', real(1/1j), imag(1/1j));
fprintf('j*z = %g %+gj  : z turned by +90 deg, same length %.1f\n\n', real(1j*z), imag(1j*z), abs(1j*z));

%% 4. atan vs atan2
zz = -3 - 4j;
fprintf('z = -3-4j: atan(b/a) = %.2f deg (wrong quadrant), atan2(b,a) = %.2f deg, angle() = %.2f deg\n', ...
        atan(imag(zz)/real(zz))*180/pi, atan2(imag(zz), real(zz))*180/pi, angle(zz)*180/pi);

%% 5. Figures: the same algebra drawn as arrows
hd = @(p0, p1, s, sg) p1 - s*(p1 - p0)/abs(p1 - p0)*exp(sg*0.45j);   % arrowhead corner
arrow = @(p0, p1, col) plot(real([p0 p1 NaN hd(p0,p1,0.45,1) p1 hd(p0,p1,0.45,-1)]), ...
                            imag([p0 p1 NaN hd(p0,p1,0.45,1) p1 hd(p0,p1,0.45,-1)]), 'Color', col, 'LineWidth', 2);
uc = exp(1j*linspace(0, 2*pi, 200));
figure('Position', [80 80 1100 360]);
subplot(1,3,1); hold on; grid on;
arrow(0, z, [0.12 0.37 0.82]); arrow(z, z+w, [0.76 0.21 0.44]); arrow(0, z+w, [0.06 0.48 0.39]);
axis equal; axis([-1 8 -1 8]); xlabel('real part'); ylabel('imaginary part');
title('z + w: tip to tail  (4 + 6j)');
subplot(1,3,2); hold on; grid on;
plot(real(uc), imag(uc), ':', 'Color', [0.5 0.5 0.5]);
arrow(0, z, [0.12 0.37 0.82]); arrow(0, w, [0.76 0.21 0.44]); arrow(0, z*w, [0.06 0.48 0.39]);
axis equal; axis([-7 6 -1 11]); xlabel('real part'); ylabel('imaginary part');
title('z w: lengths multiply, angles add (-5 + 10j)');
subplot(1,3,3); hold on; grid on;
arrow(0, z, [0.12 0.37 0.82]); arrow(0, conj(z), [0.66 0.43 0.02]); arrow(0, z/w, [0.06 0.48 0.39]);
plot([-1 6], [0 0], 'k-');
axis equal; axis([-1 6 -5 5]); xlabel('real part'); ylabel('imaginary part');
title('z* = mirror (3 - 4j);  z / w = 2.2 - 0.4j');
