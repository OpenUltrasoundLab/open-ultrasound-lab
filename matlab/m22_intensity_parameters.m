%% m22_intensity_parameters.m  -  From a hydrophone waveform to PII, t_d, I_pa, I_ta and duty factor
% Purpose : compute the exposure quantities exactly as defined for an output report:
%           p_c, p_r, PII = int p^2/(rho c) dt, pulse duration t_d = 1.25 (t90 - t10) of the
%           cumulative p^2 integral, I_pa = PII/t_d, acoustic working frequency (zero crossings
%           and -3 dB spectral mean), then I_ta and DF for scanned, PW, plane-wave and push modes.
% Change  : f0, ncyc (cycles at -6 dB of the Gaussian envelope), p0, h2 (2nd-harmonic level),
%           frame rate, beam width W6 and line spacing dx, prf.
% Expect  : 3 MHz, 3 cycles, 2 MPa (linear): PII = 0.100 mJ/cm^2; t_d = 0.905 us from the waveform
%           (0.962 us from the envelope-only closed form: the carrier ripple in p^2 shifts t10 and
%           t90 for few-cycle pulses); PW at 5 kHz 502 mW/cm^2; B-mode 30 Hz 11.3 mW/cm^2.
rng(1);
rhoc = 1.5e6;                      % Rayl, water
f0 = 3e6; ncyc = 3; p0 = 2e6;      % Hz, cycles at -6 dB, Pa
fs = 500e6; t = (-12e-6:1/fs:12e-6); % s
sig = ncyc/(2.3548*f0);            % s, envelope standard deviation
env = exp(-t.^2/(2*sig^2));
%% 1. Linear pulse and a pulse with a 2nd harmonic (nonlinear-looking, p_c > p_r)
h2 = 0.35;                         % relative 2nd-harmonic amplitude
P = [p0*env.*cos(2*pi*f0*t); p0*env.*(cos(2*pi*f0*t) + h2*env.*cos(2*pi*2*f0*t))];
names = {'linear', 'with 2nd harmonic'};
for k = 1:2
  p = P(k,:);
  pc = max(p); pr = -min(p);
  cum = cumsum(p.^2)/fs/rhoc;  PII = cum(end);              % J/m^2
  t10 = t(find(cum >= 0.1*PII, 1)); t90 = t(find(cum >= 0.9*PII, 1));
  td = 1.25*(t90 - t10);  Ipa = PII/td;
  % acoustic working frequency: zero crossings within the -12 dB envelope, and -3 dB spectral mean
  sel = find(env > 0.25); ps = p(sel); ts = t(sel);
  ic = find(ps(1:end-1).*ps(2:end) < 0);                    % sign changes
  tc = ts(ic) - ps(ic).*(ts(ic+1) - ts(ic))./(ps(ic+1) - ps(ic));  % interpolated crossing times
  fzc = (numel(tc) - 1)/(2*(tc(end) - tc(1)));
  N = 2^nextpow2(numel(p))*4; S = abs(fft(p, N)); fr = (0:N-1)*fs/N; half = 1:N/2;
  S = S(half); fr = fr(half); [Sm, im] = max(S);
  il = find(S(1:im) < Sm/sqrt(2), 1, 'last') + 1; iu = im - 1 + find(S(im:end) < Sm/sqrt(2), 1);
  fawf = (fr(il) + fr(iu))/2;
  fprintf('%s pulse:\n', names{k});
  fprintf('  p_c = %.2f MPa, p_r = %.2f MPa\n', pc/1e6, pr/1e6);
  fprintf('  PII = %.4f mJ/cm^2, t_d = %.3f us, I_pa = %.1f W/cm^2\n', PII*1e3/1e4, td*1e6, Ipa/1e4);
  fprintf('  f_awf: zero-crossing %.2f MHz, -3 dB spectral mean %.2f MHz\n', fzc/1e6, fawf/1e6);
  if k == 1
    fprintf('  closed forms: t_d = 2.2655 sigma = %.3f us, PII = p^2/(2 rho c) sigma sqrt(pi) = %.4f mJ/cm^2, I_pa = 0.782 p^2/(2 rho c) = %.1f W/cm^2\n', ...
      2.2655*sig*1e6, p0^2/(2*rhoc)*sig*sqrt(pi)*1e3/1e4, 0.7824*p0^2/(2*rhoc)/1e4);
    PIIlin = PII; tdlin = td; cumlin = cum; t10l = t10; t90l = t90;
  end
end
% envelope closed form vs exact t_d as the pulse gets longer
fprintf('\n  cycles | t_d exact / closed form\n');
for nc = [2 3 5 10 20]
  sg = nc/(2.3548*f0); q = exp(-t.^2/(2*sg^2)).*cos(2*pi*f0*t); cq = cumsum(q.^2); cq = cq/cq(end);
  tdq = 1.25*(t(find(cq >= 0.9, 1)) - t(find(cq >= 0.1, 1)));
  fprintf('   %4d  |  %.3f\n', nc, tdq/(2.2655*sg));
end
%% 2. Pulses per second at the point, I_ta and duty factor
frame = 30; W6 = 1.5e-3; dx = 0.3e-3;                 % Hz, m, m
w = W6/1.6651;                                         % 1/e^2 intensity half-width
lines = (-60:60)*dx; neff = sum(exp(-2*lines.^2/w^2)); % exact sum over lines
fprintf('\nScanned B-mode: overlap sum %.3f pulses/frame (approx 0.7527 W6/dx = %.3f)\n', neff, 0.7527*W6/dx);
modes = {'B-mode 30 Hz', 'colour 10 Hz x ens 10', 'PW 5 kHz', 'plane wave 5 kHz'};
rate = [frame*neff, 10*neff*10, 5000, 5000];
fprintf('  mode                   | pulses/s at point | I_ta (mW/cm^2) | DF\n');
for k = 1:numel(modes)
  Ita = PIIlin*rate(k);
  fprintf('  %-22s | %10.1f        | %10.2f     | %.2e\n', modes{k}, rate(k), Ita*1e3/1e4, tdlin*rate(k));
end
% push pulse: 1200 cycles at -6 dB (400 us), 3 MPa, 3 MHz, 2 per second
sp = 1200/(2.3548*3e6); PIIp = (3e6)^2/(2*rhoc)*sp*sqrt(pi);
fprintf('  push 400 us, 3 MPa, 2/s | PII = %.3f J/cm^2 | I_ta = %.0f mW/cm^2 | DF = %.2e\n', PIIp/1e4, 2*PIIp*1e3/1e4, 2*2.2655*sp);
fprintf('  push / 2-cycle PII ratio at equal pressure = %.0f\n', 1200/2);
fprintf('Track 3: I_SPTA.3 <= 720 mW/cm^2 and (MI <= 1.9 or I_SPPA.3 <= 190 W/cm^2)\n');
%% Figure
figure('Position', [100 100 900 330]);
subplot(1,2,1);
plot(t*1e6, P(1,:)/1e6, 'b', t*1e6, P(2,:)/1e6, 'm'); grid on; xlim([-1.5 1.5]);
xlabel('t (\mus)'); ylabel('p (MPa)'); legend(names); title('waveforms');
subplot(1,2,2);
plot(t*1e6, cumlin/PIIlin, 'b', 'LineWidth', 1.5); hold on;
plot([t10l t10l]*1e6, [0 1], 'k--', [t90l t90l]*1e6, [0 1], 'k--'); grid on; xlim([-1.5 1.5]);
xlabel('t (\mus)'); ylabel('cumulative \int p^2 (norm.)');
title(sprintf('t_d = 1.25 (t_{90} - t_{10}) = %.3f \\mus', tdlin*1e6));
