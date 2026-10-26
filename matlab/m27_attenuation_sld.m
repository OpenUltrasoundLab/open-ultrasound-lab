%% m27_attenuation_sld.m  -  attenuation by the spectral log difference, with a reference phantom
% Purpose : simulate speckle RF A-lines in which EVERY scatterer's echo carries
%           its own two-way attenuation exp(-2 alpha f z) and a toy focusing
%           (diffraction) factor D(f,z), for a sample (alpha_s) and a reference
%           phantom (alpha_r) scanned with the same system. Then
%           (1) show that ln S_s(f,z) - ln S_r(f,z) falls with depth at slope
%               -4 (alpha_s - alpha_r) f  (POWER spectra, two-way, Np/cm);
%           (2) estimate alpha by the reference phantom method (regression over
%               depth, then over frequency) and by the two-gate spectral log
%               difference (SLD), and show the diffraction bias you get without
%               the reference;
%           (3) measure the estimator spread over independent groups of A-lines
%               and compare with sigma_alpha ~ sqrt(1.5 g_w psi1(K)/(B^3 T))/dz
%               (Np/cm/MHz; g_w = 1.75 for a Hann gate, 0.93 rectangular,
%               Monte Carlo on white noise; psi1 = trigamma, about 1/K);
%           (4) run the centroid-downshift estimator on the same data, raw
%               and with the reference phantom's centroid slope subtracted.
% Change  : alpha_s, gate length Lg, lines per estimate K, band, focus zF, snr_dB.
% Expect  : ratio of measured to -4 x dalpha slope = 1.00; reference method
%           0.70 dB/cm/MHz; without the reference the focusing factor biases
%           the estimate to about 0.96; two-gate SLD std 0.067 dB/cm/MHz for
%           K = 32, Lg = 4 mm, dz = 2 cm, B = 4 MHz (formula 0.069); regression
%           over 31 gates (1.5-4.5 cm) std 0.030; raw centroid 1.13 (focusing
%           bias), reference-subtracted centroid 0.75 (std 0.030, 6 % model bias).
clear; close all; rng(1);

%% Parameters
c   = 0.154;            % cm/us (1540 m/s)
fs  = 40;               % MHz
f0  = 5;  bw = 0.6;     % MHz, -6 dB fractional bandwidth of the AMPLITUDE spectrum
sigf = bw*f0/(2*sqrt(2*log(2)));   % amplitude-spectrum std (MHz); power std = sigf/sqrt(2)
NP  = 20/log(10);       % 8.686 dB per neper
alpha_s = 0.70;         % sample attenuation, dB/cm/MHz
alpha_r = 0.50;         % reference phantom attenuation, dB/cm/MHz (known)
zF  = 2.0;  Fnum = 2;   % focal depth (cm) and F-number of the toy focusing model
snr_dB = Inf;           % electronic SNR at the shallowest gate (Inf = no noise)
zmax = 5.0;             % deepest scatterer (cm)
Nlines = 4096;          % independent A-lines per medium
K   = 32;               % A-lines averaged per estimate (one ROI)
Lg  = 0.4;              % gate length (cm) = 4 mm
band = [3 7];           % analysis band (MHz)

%% Simulation: RF(f) = G(f) sum_n s_n sqrt(D(f,z_n)) exp(-2 alpha f z_n) exp(-j 2 pi f t_n)
Nt = 4096;  tn = (0:ceil(2*zmax/c*fs))'/fs;  zn = c*tn/2;      % scatterer times (us) and depths (cm)
k1 = round(0.5/fs*Nt);  k2 = round(10.5/fs*Nt);  fk = (k1:k2)'*fs/Nt;   % 0.5-10.5 MHz (G < -68 dB outside)
G  = exp(-(fk - f0).^2/(2*sigf^2));                            % amplitude transfer function
zR = 3.5*(c./fk)*Fnum^2;                                       % toy half depth of field 3.5 lambda F#^2 (cm), larger at low f
Dsq = 1./sqrt(1 + ((repmat(zn', numel(fk), 1) - zF)./repmat(zR, 1, numel(zn))).^2);  % sqrt of intensity factor
Ph  = exp(-1i*2*pi*fk*tn');                                    % propagation phase
tic;
Ms = repmat(G, 1, numel(zn)).*Dsq.*exp(-2*(alpha_s/NP)*fk*zn').*Ph;   % frequency x scatterer transfer, sample
Mr = repmat(G, 1, numel(zn)).*Dsq.*exp(-2*(alpha_r/NP)*fk*zn').*Ph;   % same for the reference phantom
makeRF = @(M, S) real(ifft([zeros(k1, size(S, 2)); M*S; ...
           zeros(Nt - k2 - 1, size(S, 2))]))*2;                % one-sided spectrum -> real RF
RFs = zeros(Nt, Nlines);  RFr = RFs;
for b = 1:Nlines/256
  idx = (b-1)*256 + (1:256);
  RFs(:, idx) = makeRF(Ms, randn(numel(zn), 256));
  RFr(:, idx) = makeRF(Mr, randn(numel(zn), 256));
end
if isfinite(snr_dB)
  zq = round(2*1.0/c*fs);  p0 = mean(mean(RFs(zq:zq+100, :).^2));
  RFs = RFs + sqrt(p0*10^(-snr_dB/10))*randn(size(RFs));
  RFr = RFr + sqrt(p0*10^(-snr_dB/10))*randn(size(RFr));
end
fprintf('Simulated %d + %d A-lines of %d scatterers in %.1f s\n', Nlines, Nlines, numel(zn), toc);

%% Gated power spectra (Hann gates, zero padded)
ng = round(2*Lg/c*fs);  T = ng/fs;  Np = 512;
w  = 0.5 - 0.5*cos(2*pi*((0:ng-1)' + 0.5)/ng);
fP = (0:Np/2)'*fs/Np;  kb = fP >= band(1) & fP <= band(2);  fb = fP(kb);  B = band(2) - band(1);
zc = (1.5:0.1:4.5)';                                           % gate centres (cm)
gateS = @(RF, z, cols) mean(abs(fft(RF(round(2*z/c*fs) - ng/2 + (1:ng), cols).*repmat(w, 1, numel(cols)), Np)).^2, 2);
Ss = zeros(Np, numel(zc));  Sr = Ss;
for i = 1:numel(zc)
  Ss(:, i) = gateS(RFs, zc(i), 1:Nlines);  Sr(:, i) = gateS(RFr, zc(i), 1:Nlines);
end
fprintf('Gate %.1f mm = %d samples = %.2f us; band %g-%g MHz; Hann; %d-point FFT\n', Lg*10, ng, T, band, Np);

%% 1. The factor 4: slope of the log ratio vs depth at each frequency (all lines, ROI 1.5-4.5 cm)
sel = zc >= 1.5 & zc <= 4.5;
RS  = log(Ss(kb, sel)./Sr(kb, sel));                            % ln power ratio, Np units
Zm  = [ones(sum(sel), 1) zc(sel)];
cf  = Zm \ RS';  slope_z = cf(2, :)';                           % d ln(Ss/Sr)/dz  per frequency (1/cm)
dA  = (alpha_s - alpha_r)/NP;                                   % true Delta alpha, Np/cm/MHz
fit1 = [ones(numel(fb), 1) fb] \ slope_z;
fprintf('\n1) d/dz ln(Ss/Sr) vs f: slope %.4f /cm/MHz; -4*dalpha = %.4f  -> ratio %.3f (expect 1)\n', ...
        fit1(2), -4*dA, fit1(2)/(-4*dA));
fprintf('   (power spectra, two-way: factor 4 = 2 (round trip) x 2 (power); amplitude spectra give 2)\n');

%% 2. Estimators on all 1024 lines
aRef = alpha_r - NP*fit1(2)/4;                                  % reference phantom method
RSn  = log(Ss(kb, sel));  cfn = Zm \ RSn';
fitn = [ones(numel(fb), 1) fb] \ cfn(2, :)';  aNoRef = -NP*fitn(2)/4;
i1 = find(abs(zc - 2) < 1e-9);  i2 = find(abs(zc - 4) < 1e-9);  dz = zc(i2) - zc(i1);
y2  = log(Ss(kb, i1)./Sr(kb, i1)) - log(Ss(kb, i2)./Sr(kb, i2));
p2  = [ones(numel(fb), 1) fb] \ y2;  aSLD = alpha_r + NP*p2(2)/(4*dz);
fprintf('\n2) alpha_s true %.3f dB/cm/MHz (reference %.2f)\n', alpha_s, alpha_r);
fprintf('   reference phantom method (regression 1.5-4.5 cm): %.3f dB/cm/MHz\n', aRef);
fprintf('   two-gate SLD, z = 2 and 4 cm, with reference : %.3f dB/cm/MHz\n', aSLD);
fprintf('   same regression WITHOUT reference            : %.3f dB/cm/MHz  <- diffraction bias\n', aNoRef);

%% 3. Variance: independent groups of K lines (reference averaged over all lines)
ng_ = Nlines/K;  aG = zeros(ng_, 1);  aR = aG;  aC = aG;
kc = fP > 0.5 & fP < 12;  sigp2 = sigf^2/2;                     % centroid band; power-spectrum variance (MHz^2)
fcr = (fP(kc)'*Sr(kc, :))./sum(Sr(kc, :));  pcr = Zm \ fcr(sel)'; % reference centroid slope vs depth (all lines)
SrA1 = Sr(kb, i1);  SrA2 = Sr(kb, i2);
zs = zc(sel);  j1 = find(abs(zs - 2) < 1e-9);  j2 = find(abs(zs - 4) < 1e-9);
for g = 1:ng_
  cols = (g-1)*K + (1:K);
  Sg = zeros(Np, numel(zs));
  for i = 1:numel(zs), Sg(:, i) = gateS(RFs, zs(i), cols); end   % each gate spectrum computed once
  y  = log(Sg(kb, j1)./SrA1) - log(Sg(kb, j2)./SrA2);           % two-gate SLD
  pg = [ones(numel(fb), 1) fb] \ y;  aG(g) = alpha_r + NP*pg(2)/(4*dz);
  cg = Zm \ log(Sg(kb, :)./Sr(kb, sel))';                       % reference method, all gates
  pr = [ones(numel(fb), 1) fb] \ cg(2, :)';  aR(g) = alpha_r - NP*pr(2)/4;
  Cg = (fP(kc)'*Sg(kc, :))./sum(Sg(kc, :));                      % centroid of every gate
  pcg = Zm \ Cg';  aC(g) = alpha_r - NP*(pcg(2) - pcr(2))/(4*sigp2);
end
psi1 = sum(1./(K:1e6).^2);                                      % trigamma psi1(K) without a toolbox
gw   = 1.75;                                                    % Hann-gate correlation factor (white-noise Monte Carlo)
sTh  = NP*sqrt(1.5*gw*psi1/(B^3*T))/dz;
fprintf('\n3) %d groups of K = %d lines:\n', ng_, K);
fprintf('   two-gate SLD : mean %.3f, std %.3f dB/cm/MHz;  formula sqrt(1.5 g_w psi1(K)/(B^3 T))/dz = %.3f\n', mean(aG), std(aG), sTh);
fprintf('   regression over %d overlapping gates, 1.5-4.5 cm: mean %.3f, std %.3f dB/cm/MHz\n', sum(sel), mean(aR), std(aR));
fprintf('   centroid, reference slope subtracted: mean %.3f, std %.3f dB/cm/MHz\n', mean(aC), std(aC));

%% 4. Centroid downshift on the same sample data (no reference)
fcen = (fP(kc)'*Ss(kc, :))./sum(Ss(kc, :));
pc = Zm \ fcen(sel)';  aCen = -NP*pc(2)/(4*sigp2);
aCenRef = alpha_r - NP*(pc(2) - pcr(2))/(4*sigp2);
fprintf('\n4) centroid slope %.4f MHz/cm; alpha = -slope/(4 sigma_p^2) = %.3f dB/cm/MHz (true %.2f)\n', pc(2), aCen, alpha_s);
fprintf('   expected slope without diffraction: -4 alpha sigma_p^2 = %.4f MHz/cm\n', -4*alpha_s/NP*sigp2);
fprintf('   with the reference centroid slope subtracted: alpha = %.3f dB/cm/MHz\n', aCenRef);

%% Figures
figure('Position', [100 100 1000 380]);
subplot(1, 3, 1);
zsel = zc(sel);  ic = find(ismember(round(10*zsel), [20 25 30 35 40]));
plot(fb, RS(:, ic), 'LineWidth', 1.2); grid on;
xlabel('frequency (MHz)'); ylabel('ln(S_s/S_r)  (Np)');
legend('z = 2.0 cm', '2.5', '3.0', '3.5', '4.0', 'Location', 'southwest');
title('Log spectral ratio vs depth');
subplot(1, 3, 2);
plot(fb, slope_z, 'b.', fb, fit1(1) + fit1(2)*fb, 'r-', fb, -4*dA*fb, 'k--', 'LineWidth', 1.2); grid on;
xlabel('frequency (MHz)'); ylabel('d ln(S_s/S_r)/dz  (1/cm)');
legend('per-frequency slope', 'linear fit', '-4\Delta\alpha f (truth)', 'Location', 'southwest');
title('Slope vs f gives -4\Delta\alpha');
subplot(1, 3, 3);
plot(zc, fcen, 'b-', zc, f0 - 4*alpha_s/NP*sigp2*zc, 'k--', 'LineWidth', 1.2); grid on;
xlabel('gate depth z (cm)'); ylabel('centroid (MHz)');
legend('measured (with focusing)', 'f_0 - 4\alpha\sigma_p^2 z', 'Location', 'southwest');
title('Centroid downshift');
