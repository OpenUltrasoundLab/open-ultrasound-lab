%% m22_bioheat_toy.m  -  Is TIS = 1 about 1 degC?  Pennes steady state, duty cycle and thermal dose
% Purpose : (1) integrate the steady-state Pennes Green's function G(r) = exp(-r/L)/(4 pi k r) over
%           the heat source q = mu_a I of a collimated Gaussian beam carrying W = W_deg = 210/f mW
%           in tissue attenuating 0.3 dB/cm/MHz (all absorbed), and report the peak rise versus
%           beam radius and frequency; (2) a first-order duty-cycle heating model and CEM43.
% Change  : k (W/m/K), L (perfusion length, m), radii a, frequencies f, tau, duty, period.
% Expect  : peak rise of order 0.5-1.6 degC (within the AIUM "factor of 2" of TI = 1), weakly
%           dependent on f (absorbed power per length ~ f W = 210 mW MHz for all f), larger for
%           narrow beams. Duty cycle: period << tau follows the mean DF*dTss.
% Note    : infinite homogeneous medium, no surface boundary or transducer self-heating.
rng(1);
k = 0.5; L = 0.01;                          % W/(m K), m
radii = [0.5 1 2 4]*1e-3;                   % m, 1/e^2 intensity radius
freqs = [1 2 3.5 5 7.5];                    % MHz
dT = zeros(numel(radii), numel(freqs));
for ia = 1:numel(radii)
  a = radii(ia);
  dr = a/24; r = ((1:6*24) - 0.5)*dr;       % m
  for jf = 1:numel(freqs)
    f = freqs(jf); W = 210/f*1e-3;          % W
    mu = 0.3*log(10)/10*f*100;              % 1/m, intensity attenuation (all absorbed)
    dz = min(a/4, 5e-4); zmax = min(0.15, 8/mu + 0.05); z = ((1:round(zmax/dz)) - 0.5)*dz;
    [R, Z] = meshgrid(r, z);
    q = mu*2*W/(pi*a^2)*exp(-2*R.^2/a^2).*exp(-mu*Z);      % W/m^3
    best = 0;
    for z0 = 0:0.0025:0.04                 % search the on-axis maximum
      D = sqrt(R.^2 + (Z - z0).^2);
      v = sum(sum(q.*exp(-D/L)./(4*pi*k*D).*2*pi.*R))*dr*dz;
      best = max(best, v);
    end
    dT(ia, jf) = best;
  end
end
fprintf('Steady Pennes rise (degC) for W = 210/f mW (TIS = 1 if all power counted), k = %.1f, L = %.0f mm\n', k, L*1e3);
fprintf('  radius (mm) |'); fprintf('  %4.1f MHz', freqs); fprintf('\n');
for ia = 1:numel(radii)
  fprintf('     %4.1f     |', radii(ia)*1e3); fprintf('   %6.2f ', dT(ia,:)); fprintf('\n');
end
fprintf('  range %.2f to %.2f degC; variation from 1 to 7.5 MHz at 1 mm: %.0f %%\n', min(dT(:)), max(dT(:)), 100*(1 - dT(2,end)/dT(2,1)));
%% Convergence check (one case, finer grid)
a = 1e-3; f = 3.5; W = 210/f*1e-3; mu = 0.3*log(10)/10*f*100;
dr = a/48; r = ((1:6*48) - 0.5)*dr; dz = a/8; z = ((1:round(0.1/dz)) - 0.5)*dz; [R, Z] = meshgrid(r, z);
q = mu*2*W/(pi*a^2)*exp(-2*R.^2/a^2).*exp(-mu*Z); best = 0;
for z0 = 0:0.0025:0.04
  D = sqrt(R.^2 + (Z - z0).^2); best = max(best, sum(sum(q.*exp(-D/L)./(4*pi*k*D).*2*pi.*R))*dr*dz);
end
fprintf('  convergence (a = 1 mm, 3.5 MHz): coarse %.3f, fine %.3f degC\n', dT(2,3), best);
%% Time constants
kappa = k/(1050*3600);                      % m^2/s
fprintf('\nkappa = %.3f mm^2/s; tau = a^2/(4 kappa): a = 0.5 mm %.2f s, 1 mm %.2f s, 5 mm %.1f s; L^2/kappa = %.0f s\n', ...
  kappa*1e6, (0.5e-3)^2/(4*kappa), (1e-3)^2/(4*kappa), (5e-3)^2/(4*kappa), L^2/kappa);
%% Duty-cycle model: dT' = (dTss*s(t) - dT)/tau
dTss = 3; tau = 20; D = 0.5; Tdur = 600; dt = 0.01;
periods = [2 60];
figure('Position', [100 100 900 330]);
subplot(1,2,2); hold on;
for ip = 1:numel(periods)
  Tp = periods(ip); t = 0:dt:Tdur; s = mod(t, Tp) < D*Tp;
  th = zeros(size(t)); e = exp(-dt/tau);
  for n = 2:numel(t), th(n) = dTss*s(n) + (th(n-1) - dTss*s(n))*e; end
  T = 37 + th; R = 0.25*(T < 43) + 0.5*(T >= 43); cem = sum(dt/60*R.^(43 - T));
  fprintf('period %5.1f s: peak dT %.2f degC, mean (last half) %.2f, DF*dTss = %.2f, CEM43 = %.3g min\n', ...
    Tp, max(th), mean(th(t > Tdur/2)), D*dTss, cem);
  plot(t, th, 'LineWidth', 1.2);
end
xlabel('t (s)'); ylabel('\DeltaT (\circC)'); grid on; legend('period 2 s', 'period 60 s'); title(sprintf('\\tau = %g s, DF = %g, \\DeltaT_{ss} = %g', tau, D, dTss));
fprintf('CEM43 checks: 40 degC for 128 min = %.2f; 41 degC for 30 min = %.3f; 45 degC for 1 min = %.1f\n', 128*0.25^3, 30*0.25^2, 0.5^-2);
subplot(1,2,1);
plot(radii*1e3, dT, '-o', 'LineWidth', 1.3); hold on; plot([0 4.5], [1 1], 'k--', [0 4.5], [2 2], 'k:', [0 4.5], [0.5 0.5], 'k:');
xlabel('beam radius a (mm)'); ylabel('peak \DeltaT (\circC)'); grid on;
legend('1 MHz', '2 MHz', '3.5 MHz', '5 MHz', '7.5 MHz'); title('Pennes steady rise at W = 210/f mW');
