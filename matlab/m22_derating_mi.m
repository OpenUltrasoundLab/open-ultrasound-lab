%% m22_derating_mi.m  -  Derating at 0.3 dB/cm/MHz, the depth z_MI, and the mechanical index
% Purpose : (1) print the derating factors for pressure (20 log) and intensity (10 log);
%           (2) compute the on-axis CW field of a focused circular piston in water by a
%               numerical Rayleigh ring sum, derate it, find the depth of the derated maximum
%               (z_MI) and the MI; (3) show how MI and z_MI move with frequency for the same
%               water peak pressure, and the water pressure that MI = 1.9 allows.
% Change  : f0 (MHz), a (mm), d (mm), p_water (MPa).
% Expect  : 3 MHz at 5 cm: 4.50 dB, p x 0.596, I x 0.355. For a = 6 mm, d = 70 mm, 3.5 MHz,
%           4 MPa water peak: water peak at 46.0 mm (focal shift), z_MI = 41.5 mm, MI = 1.26.
% Notes   : Linear CW model (one frequency), so p_r^2 is proportional to PII at every depth and
%           the IEC depth (max derated PII) equals the depth of the max derated p_r.
rng(1);
c = 1500;                 % m/s, water
%% 1. Derating table
fs = [2 3 3.5 5 7.5 10];  % MHz
zs = [2 5 10];            % cm
fprintf('Derating 0.3 dB/cm/MHz:  loss dB | pressure factor | intensity factor\n');
for f = fs
  for z = zs
    A = 0.3*f*z;
    fprintf('  f = %4.1f MHz, z = %4.1f cm: %5.2f dB | p x %.3f | I x %.3f\n', f, z, A, 10^(-A/20), 10^(-A/10));
  end
end
%% 2. Focused piston, on-axis field (ring sum of the Rayleigh integral)
f0 = 3.5;  a = 6;   d = 70;  p_water = 4.0;      % MHz, mm, mm, MPa
nr = 400; r = ((1:nr) - 0.5)/nr*a*1e-3;          % m, ring radii
dr = a*1e-3/nr;
z  = (2:0.25:150)*1e-3;                          % m, depth
onaxis = @(f) abs(sum(bsxfun(@times, exp(-1i*2*pi*f*1e6/c*(sqrt(bsxfun(@plus, z.^2, r'.^2)) - (sqrt(r'.^2 + (d*1e-3)^2) - d*1e-3)*ones(size(z)))) ./ sqrt(bsxfun(@plus, z.^2, r'.^2)), r'*dr), 1));
pw = onaxis(f0); pw = p_water*pw/max(pw);         % MPa, water
zcm = z*100;
p3 = pw.*10.^(-0.3*f0*zcm/20);                    % derated pressure
pii3 = p3.^2;                                     % proportional to derated PII
[~, iw] = max(pw); [~, im] = max(pii3);
MI = p3(im)/sqrt(f0);
fprintf('\nFocused piston a = %g mm, d = %g mm, f = %.1f MHz, water peak %.1f MPa\n', a, d, f0, p_water);
fprintf('  water peak at z = %.1f mm\n', z(iw)*1e3);
fprintf('  derated maximum (z_MI) at z = %.1f mm, loss there %.2f dB\n', z(im)*1e3, 0.3*f0*zcm(im));
fprintf('  p_r.3 = %.2f MPa, MI = %.2f\n', p3(im), MI);
fprintf('  (wrong) water peak derated at its own depth: MI = %.2f\n', pw(iw)*10^(-0.3*f0*zcm(iw)/20)/sqrt(f0));
%% 3. Same geometry and water peak at other frequencies
fl = 2:0.5:7; MIf = zeros(size(fl)); zMI = MIf; zW = MIf; pAllow = MIf;
for k = 1:numel(fl)
  q = onaxis(fl(k)); q = p_water*q/max(q);
  q3 = q.*10.^(-0.3*fl(k)*zcm/20);
  [~, j] = max(q3); [~, jw] = max(q);
  MIf(k) = q3(j)/sqrt(fl(k)); zMI(k) = z(j)*1e3; zW(k) = z(jw)*1e3;
  pAllow(k) = p_water*1.9/MIf(k);              % water peak that gives MI = 1.9 (linear scaling)
end
fprintf('\n  f (MHz) | water peak z (mm) | z_MI (mm) | MI at %.1f MPa water | water peak for MI 1.9 (MPa)\n', p_water);
for k = 1:2:numel(fl)
  fprintf('   %4.1f   |      %5.1f        |   %5.1f   |        %.2f          |        %.2f\n', fl(k), zW(k), zMI(k), MIf(k), pAllow(k));
end
fprintf('\nMI = 1.9 means p_r.3 = 1.9*sqrt(f): 3 MHz %.2f, 5 MHz %.2f, 7.5 MHz %.2f MPa\n', 1.9*sqrt(3), 1.9*sqrt(5), 1.9*sqrt(7.5));
%% Figures
figure('Position', [100 100 900 360]);
subplot(1,2,1);
plot(z*1e3, pw, 'k', z*1e3, p3, 'b', 'LineWidth', 1.5); hold on;
plot(z(iw)*1e3, pw(iw), 'ko', z(im)*1e3, p3(im), 'bo', 'MarkerFaceColor', 'b');
xlabel('depth z (mm)'); ylabel('p_r (MPa)'); legend('water', 'derated 0.3 dB/cm/MHz'); grid on; xlim([0 150]);
title(sprintf('a = %g mm, d = %g mm, %.1f MHz: MI = %.2f', a, d, f0, MI));
subplot(1,2,2);
plot(fl, MIf, 'b-o', 'LineWidth', 1.5); hold on; plot(fl, 1.9*ones(size(fl)), 'r--');
xlabel('frequency (MHz)'); ylabel('MI at fixed water peak'); grid on;
title(sprintf('water peak %.1f MPa, same geometry', p_water));
