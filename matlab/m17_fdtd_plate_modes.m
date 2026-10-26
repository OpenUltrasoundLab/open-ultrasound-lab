%% m17_fdtd_plate_modes.m  --  2D elastic FDTD of a free piezo-like bar: lateral modes vs w/t
% Purpose : show why 1D (KLM) models fail for array elements. A bar of width w
%           and thickness t (plane strain: the element is long in elevation) is
%           driven by a uniform thickness-direction eigenstress, like the field
%           between two electrodes. The spectrum of the mean face velocity
%           shows the thickness mode, the width modes, and their coupling.
% Method  : velocity-stress staggered grid (Virieux 1986), isotropic solid,
%           nu = 0.32 (ceramic-like), free surfaces by the vacuum formulation,
%           16 cells per thickness, Courant number 0.5. Normalised units:
%           t = 1, cL = 1, rho = 1, so f t / cL = 0.5 is the plate thickness
%           mode (5 MHz for a 0.456 mm plate with cL = 4560 m/s).
% Change  : wt list, nu, ppt (cells per thickness), Tend.
% Expect  : w/t = 0.2 -> single resonance at f t/cL ~ 0.44 (thin bar, theory
%           0.441); w/t ~ 0.8-2 -> split, coupled modes; w/t = 6 -> ~0.5 with
%           width harmonics. First width mode follows 0.441/(w/t).
%           Runtime ~ 15 s in Octave.
clear; close all; rng(1);

nu = 0.32; ppt = 16; Tend = 60; damp = 0.002;
cL = 1; cS = cL*sqrt((1 - 2*nu)/(2*(1 - nu))); mu = cS^2; lam = cL^2 - 2*mu;
dx = 1/ppt; dt = 0.5*dx/cL; nt = round(Tend/dt); g = 1 - damp;
sig = 0.35/0.5; t0 = 4*sig; gp = @(u) exp(-0.5*((u - t0)/sig).^2);   % Gaussian eigenstress drive
fr = 0.02:0.005:1.0;                                  % analysis frequencies (f t / cL)
win = 0.5 - 0.5*cos(2*pi*(0:nt-1)/nt);                % Hann window on the record
E = exp(-1i*2*pi*(0:nt-1)'*dt*fr);                    % DFT matrix (nt x nf)

wts = [0.2 0.3:0.1:3.0 4 6];
chart = []; specs = zeros(numel(wts), numel(fr));
tic;
for q = 1:numel(wts)
  wt = wts(q); nxb = max(2, round(wt*ppt)); nzb = ppt; pad = 2;
  nx = nxb + 2*pad; nz = nzb + 2*pad;
  S = false(nz, nx); S(pad+1:pad+nzb, pad+1:pad+nxb) = true;          % solid cells (rows = z, cols = x)
  Sx = [S(:, 2:end), false(nz, 1)]; Sz = [S(2:end, :); false(1, nx)]; Sxz = [Sz(:, 2:end), false(nz, 1)];
  rx = S + Sx; bx = zeros(nz, nx); bx(rx > 0) = 2./rx(rx > 0);       % buoyancy at vx nodes
  rz = S + Sz; bz = zeros(nz, nx); bz(rz > 0) = 2./rz(rz > 0);       % buoyancy at vz nodes
  Mxz = mu*(S & Sx & Sz & Sxz);                                      % shear modulus at sxz nodes
  L2M = (lam + 2*mu)*S; LL = lam*S;
  vx = zeros(nz, nx); vz = vx; sxx = vx; szz = vx; sxz = vx; rec = zeros(1, nt);
  I = 2:nx-1; J = 2:nz-1; top = pad + nzb;                           % row of the top face vz nodes
  for n = 1:nt
    tt = (n - 1)*dt; sv = gp(tt) - gp(tt - dt);
    vx(J,I) = g*vx(J,I) + dt/dx*bx(J,I).*(sxx(J,I+1) - sxx(J,I) + sxz(J,I) - sxz(J-1,I));
    vz(J,I) = g*vz(J,I) + dt/dx*bz(J,I).*(sxz(J,I) - sxz(J,I-1) + szz(J+1,I) - szz(J,I));
    ex = (vx(J,I) - vx(J,I-1))/dx; ez = (vz(J,I) - vz(J-1,I))/dx;
    sxx(J,I) = sxx(J,I) + dt*(L2M(J,I).*ex + LL(J,I).*ez);
    szz(J,I) = szz(J,I) + dt*(LL(J,I).*ex + L2M(J,I).*ez) + sv*S(J,I);
    sxz(J,I) = sxz(J,I) + dt*Mxz(J,I).*((vx(J+1,I) - vx(J,I)) + (vz(J,I+1) - vz(J,I)))/dx;
    rec(n) = mean(vz(top, pad+1:pad+nxb));
  end
  sp = abs((rec.*win)*E); sp = sp/max(sp); specs(q,:) = sp;
  pk = find(sp(2:end-1) > sp(1:end-2) & sp(2:end-1) >= sp(3:end) & sp(2:end-1) > 0.1) + 1;
  chart = [chart; [wt*ones(numel(pk), 1), fr(pk)', sp(pk)']]; %#ok<AGROW>
  band = pk(fr(pk) > 0.35 & fr(pk) < 0.62); [~, ib] = max(sp(band));
  if any(abs(wt - [0.2 0.6 1 1.5 2 3 6]) < 1e-9)
    fprintf('w/t = %.1f: peaks at f t/cL = %s (strongest in 0.35-0.62: %.3f)\n', wt, sprintf('%.3f ', fr(pk)), fr(band(ib)));
  end
end
fprintf('thin-bar limit (plane strain) f t/cL = 0.5*sqrt((1-2nu))/(1-nu) = %.3f; plate = 0.500\n', 0.5*sqrt(1 - 2*nu)/(1 - nu));
fprintf('first width mode estimate 0.441/(w/t): w/t = 1.5 -> %.3f, 2 -> %.3f, 3 -> %.3f\n', 0.441/1.5, 0.441/2, 0.441/3);
fprintf('total FDTD time %.1f s for %d bars\n', toc, numel(wts));

%% Figures
figure('Position', [60 60 1000 400]);
subplot(1,2,1); sel = find(ismember(round(wts*10), [2 6 10 15 30]));
for k = 1:numel(sel), plot(fr, specs(sel(k),:) - 1.1*(k-1), 'LineWidth', 1.3); hold on; end
grid on; xlabel('f t / c_L  (0.5 = 5 MHz for t = 0.456 mm)'); ylabel('normalised |v_z| (offset)');
legend(arrayfun(@(v) sprintf('w/t = %.1f', v), wts(sel), 'UniformOutput', false), 'Location', 'east');
title('Mean face velocity spectrum');
subplot(1,2,2); scatter(chart(:,1), chart(:,2), 8 + 40*chart(:,3), chart(:,3), 'filled'); hold on;
ww = linspace(0.15, 6, 300); plot(ww, 0.441./ww, 'k--', ww, 3*0.441./ww, 'k:', [0 6], [0.5 0.5], 'r-', [0 6], [0.441 0.441], 'b-');
ylim([0 1]); xlim([0 6.2]); grid on; xlabel('width / thickness  w/t'); ylabel('resonance  f t / c_L');
title('Mode chart: KLM (red line) only knows the plate'); colorbar;
