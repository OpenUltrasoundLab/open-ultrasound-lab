%% m27_bsc_esd_fit.m  -  backscatter coefficient shape, Gaussian form factor, effective scatterer diameter
% Purpose : speckle gates are simulated from a scatterer field whose power
%           spectrum is f^4 x F(f) with the Gaussian form factor
%           F = exp(-0.827 k^2 a^2) (a = effective radius, k = 2 pi f / c;
%           Insana et al. 1990), attenuated two-way to depth z, and seen through
%           a 5 MHz, 60 % system. A reference phantom with known BSC shape f^4
%           (Rayleigh scatterers) and known attenuation is scanned with the same
%           system. Then
%           (1) the reference phantom ratio gives BSC_s(f)/BSC_r(f) after
%               attenuation compensation exp(+4 (alpha_s - alpha_r) f z);
%           (2) ln(BSC_s / f^4) is linear in k^2 with slope -0.827 a^2, so a
%               linear fit returns the ESD = 2a;
%           (3) the ESD is identifiable only when ka over the band reaches ~0.5-1:
%               small scatterers give a nearly flat F and the fit fails;
%           (4) an attenuation-compensation error of 0.1 dB/cm/MHz at 3 cm
%               shifts the ESD by tens of micrometres (they share the k^2 / f slope).
% Change  : band, gate length Lg, lines per estimate K, depth z, alpha error.
% Expect  : ESD unbiased within 1 um and std 2-7 um once ka_max >= 1.4
%           (ESD >= 100 um at 7 MHz); ESD = 20 um: 36 % of fits give a slope of
%           the wrong sign and the rest average 33 um (biased, std 13 um); a
%           +0.1 dB/cm/MHz compensation error at z = 3 cm turns 100 um into 78 um.
clear; close all; rng(1);

%% Parameters
c   = 1540;                 % m/s  (= um/us, so k = 2 pi f / c in rad/um with f in MHz)
cc  = 0.154;                % cm/us, for gate lengths
fs  = 40;  f0 = 5;  bw = 0.6;
sigf = bw*f0/(2*sqrt(2*log(2)));
NP  = 20/log(10);           % 8.686 dB/Np
z   = 3.0;                  % gate depth (cm)
alpha_s = 0.5;  alpha_r = 0.5;   % dB/cm/MHz (equal here so attenuation errors are isolated)
Lg  = 0.4;  ng = round(2*Lg/cc*fs);  T = ng/fs;    % 4 mm gate
K   = 32;   Kref = 1024;    % sample lines per estimate, reference lines
band = [3 7];               % MHz
ntr = 200;                  % Monte Carlo estimates per ESD
Nf  = 1024;  ff = (0:Nf-1)'*fs/Nf;  ff(ff >= fs/2) = ff(ff >= fs/2) - fs;  fa = abs(ff);
w   = 0.5 - 0.5*cos(2*pi*((0:ng-1)' + 0.5)/ng);
Np  = 512;  fP = (0:Np/2)'*fs/Np;  kb = fP >= band(1) & fP <= band(2);  fb = fP(kb);  kk = 2*pi*fb/c;
G   = exp(-(fa - f0).^2/(2*sigf^2));
% amplitude filter for a gate: system x f^2 (Rayleigh) x sqrt(form factor) x two-way attenuation
Hamp = @(a_um, alphadB) G.*(fa/f0).^2.*exp(-0.827*(2*pi*fa/c*a_um).^2/2).*exp(-2*alphadB/NP*fa*z);
gateSpec = @(H, n) mean(abs(fft(real(ifft(fft(randn(Nf, n)).*repmat(H, 1, n)))...
                   (1:ng, :).*repmat(w, 1, n), Np)).^2, 2);

%% Reference phantom spectrum (known BSC shape f^4, known alpha_r), many lines
Sr = gateSpec(Hamp(0, alpha_r), Kref);  Sr = Sr(kb);
fprintf('Gate %.1f mm at z = %.1f cm, band %g-%g MHz, K = %d lines/estimate, reference %d lines\n', Lg*10, z, band, K, Kref);

%% 1-3. ESD sweep
ESDs = [20 40 60 80 100 150 200 300];          % um
fprintf('\n ESD(um)  ka range     mean ESD^  std    fails   ln-BSC slope n at f0\n');
X = [ones(numel(kk), 1) kk.^2];
keep = zeros(numel(fb), 3);
for e = 1:numel(ESDs)
  a = ESDs(e)/2;  est = nan(ntr, 1);
  for t = 1:ntr
    Ss = gateSpec(Hamp(a, alpha_s), K);
    R  = Ss(kb)./Sr.*exp(4*(alpha_s - alpha_r)/NP*fb*z);   % = BSC_s/BSC_r = F(f) x const
    p  = X \ log(R);                                       % ln R = b0 - 0.827 a^2 k^2
    a2 = -p(2)/0.827;
    if a2 > 0, est(t) = 2*sqrt(a2); end
    if t == 1 && any(ESDs(e) == [20 100 300]), keep(:, find(ESDs(e) == [20 100 300])) = log(R) - p(1); end
  end
  n_f0 = 4 - 2*0.827*(2*pi*f0/c*a)^2;                     % local power-law exponent of f^4 F at f0
  fprintf('  %4d    %.2f-%.2f   %7.1f  %6.1f   %5.3f   %6.2f\n', ESDs(e), kk(1)*a, kk(end)*a, ...
          mean(est(~isnan(est))), std(est(~isnan(est))), mean(isnan(est)), n_f0);
end
fprintf('BSC = f^4 exp(-0.827 k^2 a^2) peaks at ka = sqrt(2/0.827) = %.3f; local exponent n = 4 - 1.654 (ka)^2\n', sqrt(2/0.827));

%% 4. Attenuation-compensation error leaks into the ESD (noise-free spectra)
fprintf('\nCompensation error (assumed alpha_s - true), z = %.1f cm, noise-free expected spectra:\n', z);
Hn = @(a) (Hamp(a, alpha_s)).^2;                        % expected power (no window smoothing)
for a_true = [50 100]
  for dA = [-0.1 0.05 0.1 0.2]
    Pn = Hn(a_true);  Pr = Hamp(0, alpha_r).^2;
    Pn = interp1(ff(1:Nf/2), Pn(1:Nf/2), fb);  Pr = interp1(ff(1:Nf/2), Pr(1:Nf/2), fb);
    R  = Pn./Pr.*exp(4*(alpha_s + dA - alpha_r)/NP*fb*z);
    p  = X \ log(R);  a2 = -p(2)/0.827;
    if a2 > 0, s = sprintf('%.1f um', 2*sqrt(a2)); else, s = 'no real solution'; end
    fprintf('  true ESD %3d um, error %+5.2f dB/cm/MHz -> ESD^ = %s\n', 2*a_true, dA, s);
  end
end

%% Figures
figure('Position', [100 100 950 380]);
subplot(1, 2, 1);
kaa = linspace(0, 4, 400);
plot(kaa, 10*log10(kaa.^4.*exp(-0.827*kaa.^2)/max(kaa.^4.*exp(-0.827*kaa.^2))), 'k', 'LineWidth', 1.4); hold on;
plot(kaa, 10*log10(exp(-0.827*kaa.^2)), 'b--', 'LineWidth', 1.2);
for a = [10 50 150], r = 2*pi*band/c*a; plot(r, [-29 -29], 'LineWidth', 4); end
ylim([-30 2]); grid on; xlabel('ka'); ylabel('dB');
legend('(ka)^4 F(ka), normalized', 'F = exp(-0.827 (ka)^2)', '3-7 MHz, ESD 20 um', 'ESD 100 um', 'ESD 300 um', 'Location', 'southeast');
title('Gaussian form factor and the band in ka');
subplot(1, 2, 2);
plot(kk.^2*1e3, keep, '.', 'MarkerSize', 8); hold on;
for a = [10 50 150], plot(kk.^2*1e3, -0.827*a^2*kk.^2, 'k-'); end
grid on; xlabel('k^2  (10^{-3} rad^2/um^2)'); ylabel('ln(BSC_s/BSC_r) - intercept  (Np)');
legend('ESD 20 um', 'ESD 100 um', 'ESD 300 um', 'Location', 'southwest');
title('One estimate (K lines): slope = -0.827 a^2');
