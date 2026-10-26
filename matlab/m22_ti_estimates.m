%% m22_ti_estimates.m  -  MI, derated intensities and TIS / TIB / TIC for a defined mode
% Purpose : the Gaussian-beam output model of Demo 22.4. A separable Gaussian source (1/e amplitude
%           half-width 0.4 x aperture in each axis) focused in azimuth (Fx) and elevation (Fy)
%           radiates a Gaussian-envelope pulse in water. From it: p_r.3(z), PII.3(z), I_ta.3(z),
%           I_pa.3(z), P_alpha(z); MI at z_MI (max derated PII); I_SPTA.3, I_SPPA.3; TIS (scanning:
%           P1 f/210; non-scanning large aperture: max over z >= z_bp of min[P_alpha, I_ta,alpha x 1 cm^2]
%           f/210), TIB (model P_alpha/(40 mW/cm x d_-6), scanning TIB = TIS), TIC = P/(40 D_eq).
% Change  : the preset structs (f MHz, Lx Ly mm, Fx Fy mm, p0 MPa source pressure, N cycles at -6 dB,
%           mode 'scan' | 'unscan' | 'plane', rate or lines/frame/fov).
% Expect  : PW preset: MI 0.60 at 4.80 cm, I_SPTA.3 428 mW/cm^2, P 39.2 mW, TIS 0.34 (z_s 2.36 cm),
%           TIB 1.49, TIC 0.63. B-mode: MI 1.16, TIS = TIB 0.12, TIC 0.32. Push: MI 1.62,
%           I_SPTA.3 159 mW/cm^2, I_SPPA.3 206 W/cm^2. sqrt-form TIB cross-check (C_TIB = 50):
%           1.063 x the diameter form for these Gaussian beams.
% Note    : TIB here uses the diameter model (40 mW/cm, Nowicki 2020 eq. 6). IEC 62359 computes
%           non-scanning TIB as sqrt(P_alpha I_ta,alpha)/C_TIB at the depth where the product is
%           largest. C_TIB = 50 mW/cm is the commonly quoted ODS value but was NOT confirmed from a
%           public source for this course; the script prints that form only as a cross-check.
%           Use the standard's own text for submission numbers.
rhoc = 1.5e6; cw = 1500;
zs = 0.1:0.02:15;                                     % cm
der = @(f, z, pw) 10.^(-0.3*f*z/(10*pw + 20*(1-pw)));  % pw = 1: intensity (10 log), 0: pressure
den = @(ws, F, lam) (1/F)^2 + (lam/(pi*ws^2))^2;              % complex beam parameter q = q0 + z
qr0 = @(ws, F, lam) (-1/F)/den(ws, F, lam);                     % Re q0 (converging: R = -F)
qi0 = @(ws, F, lam) (lam/(pi*ws^2))/den(ws, F, lam);            % Im q0
gw  = @(ws, F, lam, z) sqrt(lam*((qr0(ws, F, lam) + z).^2 + qi0(ws, F, lam)^2)/(pi*qi0(ws, F, lam)));  % 1/e half-width (m)
% presets (same as the page)
P1 = struct('name','PW Doppler','f',3.0,'Lx',19.2,'Ly',10,'Fx',50,'Fy',60,'p0',0.45,'N',8,'mode','unscan','rate',6000,'lines',0,'frame',0,'fov',0);
P2 = struct('name','B-mode','f',3.5,'Lx',19.2,'Ly',10,'Fx',50,'Fy',60,'p0',0.9,'N',2,'mode','scan','rate',0,'lines',128,'frame',40,'fov',38.4);
P3 = struct('name','push 2/s','f',3.0,'Lx',19.2,'Ly',10,'Fx',40,'Fy',60,'p0',1.0,'N',1200,'mode','unscan','rate',2,'lines',0,'frame',0,'fov',0);
P4 = struct('name','plane 5 kHz','f',5,'Lx',38.4,'Ly',5,'Fx',1e6,'Fy',20,'p0',0.8,'N',2,'mode','plane','rate',5000,'lines',0,'frame',0,'fov',0);
presets = {P1, P2, P3, P4};
out = cell(1, numel(presets));
for k = 1:numel(presets)
  P = presets{k};
  lam = cw/(P.f*1e6); wsx = 0.4*P.Lx*1e-3; wsy = 0.4*P.Ly*1e-3;
  wx = gw(wsx, P.Fx*1e-3, lam, zs*1e-2); wy = gw(wsy, P.Fy*1e-3, lam, zs*1e-2);   % m
  amp = sqrt(wsx*wsy./(wx.*wy));
  sig = P.N/(2.3548*P.f*1e6); Ee = sig*sqrt(pi); td = 2.2655*sig;
  p0 = P.p0*1e6; E = p0^2/(2*rhoc)*Ee*pi*wsx*wsy/2;                         % J per pulse
  if strcmp(P.mode, 'scan')
    pps = P.lines*P.frame; dx = P.fov*1e-3/P.lines; rate = zeros(size(zs));
    for i = 1:numel(zs)
      L = ceil(4*wx(i)/dx); l = -L:L; rate(i) = P.frame*sum(exp(-2*(l*dx).^2/wx(i)^2));
    end
  else
    pps = P.rate; rate = pps*ones(size(zs));
  end
  Pw = E*pps*1e3;                                                            % mW
  p = p0*amp; pii = p.^2/(2*rhoc)*Ee;
  pr3 = p.*der(P.f, zs, 0); pii3 = pii.*der(P.f, zs, 1); ita3 = pii.*rate.*der(P.f, zs, 1);
  ipa3 = pii/td.*der(P.f, zs, 1); P3a = Pw*der(P.f, zs, 1);
  [~, iMI] = max(pii3); MI = pr3(iMI)/1e6/sqrt(P.f);
  [ispta3, iT] = max(ita3*0.1); isppa3 = max(ipa3)/1e4;                   % mW/cm^2, W/cm^2
  tisA = P3a*P.f/210; tisB = ita3*0.1*P.f/210; tib = P3a./(40*1.6651*sqrt(wx.*wy)*100);
  Aap = P.Lx*P.Ly/100;
  if ~strcmp(P.mode, 'unscan')
    if strcmp(P.mode, 'scan'), fovcm = P.fov/10; else, fovcm = P.Lx/10; end
    Pb = Pw*min(1, 1/fovcm); TIS = Pb*P.f/210; TIB = TIS; Deq = sqrt(4*(fovcm*P.Ly/10)/pi); TIC = Pw/(40*Deq); zsv = NaN; zbv = NaN; zbp = NaN;
  else
    Deq = sqrt(4*Aap/pi); zbp = 1.5*Deq; Pb = NaN;
    if Aap <= 1
      TIS = Pw*P.f/210; zsv = NaN;
    else
      m = min(tisA, tisB); m(zs < zbp) = 0; [TIS, j] = max(m); zsv = zs(j);
    end
    [TIB, jb] = max(tib); zbv = zs(jb); TIC = Pw/(40*Deq);
    CTIB = 50;                                                   % mW/cm, UNVERIFIED (see header)
    [TIBsq, jq] = max(sqrt(P3a.*ita3*0.1)/CTIB);                 % sqrt(mW x mW/cm^2) = mW/cm
  end
  fprintf('%-12s: MI %.2f at z_MI %.2f cm | I_SPTA.3 %.1f mW/cm^2 (z %.2f) | I_SPPA.3 %.1f W/cm^2 | P %.1f mW\n', ...
    P.name, MI, zs(iMI), ispta3, zs(iT), isppa3, Pw);
  fprintf('              TIS %.3f (z_s %.2f, z_bp %.2f) | TIB %.3f (z_b %.2f) | TIC %.3f (D_eq %.2f cm)\n', TIS, zsv, zbp, TIB, zbv, TIC, Deq);
  if strcmp(P.mode, 'unscan')
    fprintf('              cross-check sqrt(P_a I_ta,a)/50: TIB %.3f at %.2f cm (ratio to diameter form %.3f)\n', TIBsq, zs(jq), TIBsq/TIB);
  end
  ok = ispta3 <= 720 && (MI <= 1.9 || isppa3 <= 190);
  if ok, fprintf('              Track 3: within\n'); else, fprintf('              Track 3: exceeded\n'); end
  out{k} = struct('pr3', pr3, 'tisA', tisA, 'tisB', tisB, 'tib', tib, 'zbp', zbp);
end
%% TIB and MI versus azimuth focal depth (PW preset)
Fl = [20 30 40 50 70 100];
fprintf('\nPW preset, azimuth focus sweep:\n  Fx (mm) |  MI  | TIS  | TIB\n');
for F = Fl
  P = P1; P.Fx = F;
  lam = cw/(P.f*1e6); wsx = 0.4*P.Lx*1e-3; wsy = 0.4*P.Ly*1e-3;
  wx = gw(wsx, P.Fx*1e-3, lam, zs*1e-2); wy = gw(wsy, P.Fy*1e-3, lam, zs*1e-2);
  amp = sqrt(wsx*wsy./(wx.*wy)); sig = P.N/(2.3548*P.f*1e6); Ee = sig*sqrt(pi);
  p0 = P.p0*1e6; E = p0^2/(2*rhoc)*Ee*pi*wsx*wsy/2; Pw = E*P.rate*1e3;
  pii = (p0*amp).^2/(2*rhoc)*Ee; pii3 = pii.*der(P.f, zs, 1); [~, iMI] = max(pii3);
  MI = p0*amp(iMI)*der(P.f, zs(iMI), 0)/1e6/sqrt(P.f);
  P3a = Pw*der(P.f, zs, 1); tisA = P3a*P.f/210; tisB = pii3*P.rate*0.1*P.f/210;
  Deq = sqrt(4*(P.Lx*P.Ly/100)/pi); m = min(tisA, tisB); m(zs < 1.5*Deq) = 0;
  fprintf('   %4d   | %.2f | %.2f | %.2f\n', F, MI, max(m), max(P3a./(40*1.6651*sqrt(wx.*wy)*100)));
end
%% Figure
o = out{1};
figure('Position', [100 100 900 330]);
subplot(1,2,1); plot(zs, o.pr3/1e6, 'b', 'LineWidth', 1.5); grid on; xlabel('z (cm)'); ylabel('p_{r.3} (MPa)'); title('PW preset: derated pressure');
subplot(1,2,2); plot(zs, o.tisA, 'Color', [0.8 0.5 0], 'LineWidth', 1.3); hold on; plot(zs, o.tisB, 'Color', [0 0.6 0.5], 'LineWidth', 1.3);
plot(zs, o.tib, 'm', 'LineWidth', 1.3); plot([o.zbp o.zbp], [0 1.1*max(o.tib)], 'k--'); grid on;
xlabel('z (cm)'); ylabel('index units'); legend('P_\alpha f/210', 'I_{ta,\alpha} 1 cm^2 f/210', 'TIB(z) model', 'z_{bp}'); title('TI integrands');
