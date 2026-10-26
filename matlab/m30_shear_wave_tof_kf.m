%% m30_shear_wave_tof_kf.m
% Purpose : Shear-wave speed estimators on one simulated push. Cylindrical
%           wave from a 200 us line push (sigma = 0.35 mm) in a Kelvin-Voigt
%           medium, far-field 2D Green's function
%             U(r,w) ~ S(w) exp(-j k r) / (mu* sqrt(k r)),  k = w sqrt(rho/mu*)
%           on an FFT grid (frame rate 10 kHz, 512 frames), particle velocity.
%           Estimators from 3 to 11 mm: time-to-peak (TTP) regression,
%           (peak searched within +-2 ms of the previous trace's peak),
%           cross-correlation (XC) delay regression, and the k-f phase
%           velocity (2D FFT, Hann window in x). Monte Carlo over noise.
% Change  : c0, eta, snrdB, x2 (window end), fr (frame rate), nMC.
% Expect  : eta = 0: all estimators return 1.50 m/s. eta = 0.4 Pa s:
%           c_ph(200 Hz) = 1.528 m/s, c_g(200 Hz) = 1.584 m/s; TTP on
%           particle velocity reads above c_ph (about 1.6-1.9 m/s), XC close
%           to c_ph at the trace centroid frequency (std over 40 noise draws
%           at 30 dB printed for each estimator).
rng(1);
rho = 1000; c0 = 1.5; fr = 10e3; Nt = 512; Tp = 200e-6; sr = 0.35e-3;
dx = 0.2e-3; x = (1e-3:dx:20e-3)'; nx = numel(x);
x1 = 3e-3; x2 = 11e-3; snrdB = 30; nMC = 40;
m = (1:Nt/2)'; f = m*fr/Nt; w = 2*pi*f;
etaList = [0 0.4];
for eta = etaList
  mu = rho*c0^2; ms = mu + 1j*w*eta; k = w.*sqrt(rho./ms);
  if any(imag(k) > 0), k = conj(k); end            % choose the decaying root
  S = (1 - exp(-1j*w*Tp))./(1j*w); hp = 1 - exp(-(f/30).^2); sp = exp(-real(k).^2*sr^2/2);
  V = zeros(Nt, nx);
  for ix = 1:nx
    G = exp(-1j*k*x(ix) + 1j*pi/4)./sqrt(k*x(ix))./ms;
    Z = 1j*w.*S.*G.*hp.*sp;                         % particle velocity spectrum
    full = zeros(Nt,1); full(2:Nt/2) = Z(1:end-1); full(Nt/2+1) = real(Z(end)); full(Nt/2+2:end) = conj(Z(end-1:-1:1));
    V(:,ix) = real(ifft(full));
  end
  sel = find(x >= x1 - 1e-9 & x <= x2 + 1e-9); xs = x(sel);
  cph = @(ff) 2*pi*ff./real(2*pi*ff.*sqrt(rho./(mu + 1j*2*pi*ff*eta)));
  hh = 0.5; cg = @(ff) 2*pi*2*hh./(real(2*pi*(ff+hh).*sqrt(rho./(mu + 1j*2*pi*(ff+hh)*eta))) - real(2*pi*(ff-hh).*sqrt(rho./(mu + 1j*2*pi*(ff-hh)*eta))));
  P = abs(fft(V(:, sel(1)))).^2; fcen = sum(P(2:Nt/2).*f(1:end-1))/sum(P(2:Nt/2));
  pk = max(abs(V(:, sel(1))));
  cT = zeros(nMC,1); cX = cT; cK = cT;
  for it = 1:nMC
    Vn = V(:, sel) + pk*10^(-snrdB/20)*randn(Nt, numel(sel));
    % TTP
    tt = zeros(numel(sel),1);
    for q = 1:numel(sel)
      if q == 1, rg = 1:round(0.01*fr); else, rg = max(1, imPrev - round(0.002*fr)):min(Nt, imPrev + round(0.002*fr)); end
      [~, im] = max(Vn(rg,q)); im = im + rg(1) - 1; imPrev = im; d = 0;
      if im > 1 && im < Nt, a = Vn(im-1,q); b = Vn(im,q); g = Vn(im+1,q); den = a - 2*b + g; if den < 0, d = 0.5*(a - g)/den; end, end
      tt(q) = (im - 1 + d)/fr;
    end
    pT = [xs ones(size(xs))]\tt; cT(it) = 1/pT(1);
    % XC versus first trace, FFT with zero padding
    NF = 1024; R0 = fft(Vn(:,1), NF); lg = zeros(numel(sel),1);
    for q = 1:numel(sel)
      cc = real(ifft(fft(Vn(:,q), NF).*conj(R0))); maxL = round(0.03*fr); [~, im] = max(cc(1:maxL+1)); d = 0;
      if im > 1, a = cc(im-1); b = cc(im); g = cc(im+1); den = a - 2*b + g; if den < 0, d = 0.5*(a - g)/den; end, end
      lg(q) = (im - 1 + d)/fr;
    end
    pX = [xs ones(size(xs))]\lg; cX(it) = 1/pX(1);
    % k-f at 200 Hz bin
    nW = numel(sel); win = 0.5 - 0.5*cos(2*pi*((1:nW)' - 0.5)/nW);
    TF = fft(Vn); mb = round(200/(fr/Nt)) + 1; row = TF(mb, :).'.*win;
    NK = 256; K = abs(fft(row, NK)); Kn = [K(1); K(end:-1:2)];   % index q -> k = q*dk (wave e^{-jkx})
    [~, iq] = max(Kn(1:NK/2)); d = 0; if iq > 1, a = Kn(iq-1); b = Kn(iq); g = Kn(iq+1); den = a - 2*b + g; if den < 0, d = 0.5*(a - g)/den; end, end
    kk = (iq - 1 + d)*2*pi/(NK*dx); cK(it) = 2*pi*(mb - 1)*fr/Nt/kk;
  end
  fprintf('\neta = %.1f Pa s (c0 = %.2f m/s), window %.0f-%.0f mm, %d dB, %d draws\n', eta, c0, x1*1e3, x2*1e3, snrdB, nMC);
  fprintf('  centroid frequency of the trace at 3 mm: %.0f Hz; c_ph = %.3f, c_g = %.3f m/s there\n', fcen, cph(fcen), cg(fcen));
  fprintf('  theory: c_ph(200 Hz) = %.3f m/s, c_g(200 Hz) = %.3f m/s\n', cph(200), cg(200));
  fprintf('  TTP  : %.3f +- %.3f m/s\n', mean(cT), std(cT));
  fprintf('  XC   : %.3f +- %.3f m/s\n', mean(cX), std(cX));
  fprintf('  k-f  : c_ph(200 Hz) = %.3f +- %.3f m/s\n', mean(cK), std(cK));
  if eta > 0
    figure(1); subplot(1,2,1); imagesc(x*1e3, (0:Nt-1)/fr*1e3, V.*sqrt(x')); axis xy; ylim([0 20]); xlabel('x (mm)'); ylabel('t (ms)'); title('particle velocity, \eta = 0.4 Pa s');
    subplot(1,2,2); ff = 20:5:1000; plot(ff, cph(ff), 'k', ff, cg(ff), 'r--'); xlabel('f (Hz)'); ylabel('speed (m/s)'); legend('c_{ph}', 'c_g', 'location', 'northwest'); title('Kelvin-Voigt dispersion');
  end
end
% precision scaling of (30.11)
c = 2; sig_t = 50e-6; L = 10e-3; N = 51; dxs = L/(N-1);
fprintf('\n(30.11): c = 2 m/s, sigma_t = 50 us, L = 10 mm, N = 51: sigma_c/c = %.2f %%\n', 100*c*sig_t/dxs*sqrt(12/(N*(N^2-1))));
