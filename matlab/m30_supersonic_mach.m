%% m30_supersonic_mach.m
% Purpose : Supersonic shear source. 2D finite-difference (leapfrog) solution
%           of the scalar elastic shear-wave equation in the imaging plane:
%           pushes 2 mm apart along the beam axis, each lasting 2 mm / V, so
%           the source moves at V. Measures the lateral time-to-peak speed on
%           the row crossed by the plane part of the front and compares with
%           c_s / cos(phi), sin(phi) = 1/M (30.9).
% Change  : cs, Vlist, n (number of foci).
% Expect  : c_s = 2 m/s, 5 foci: V = 15 m/s (M = 7.5, phi = 7.7 deg) reads
%           about 2.0 m/s; V = 5 m/s (M = 2.5, phi = 23.6 deg) about 2.2 m/s;
%           V = 3 m/s (M = 1.5, phi = 41.8 deg) about 2.7-2.8 m/s.
rng(1);
rho = 1000; cs = 2; n = 5; dx = 0.2e-3; nx = 121; nz = 141; x0 = -12e-3; z0 = 4e-3;
X = x0 + (0:nx-1)*dx; Zc = z0 + (0:nz-1)'*dx; [XX, ZZ] = meshgrid(X, Zc);
tRec = 0.1e-3; nRec = 80; sp = 14;
id = min(min(repmat(0:nx-1, nz, 1), repmat(nx-1:-1:0, nz, 1)), min(repmat((0:nz-1)', 1, nx), repmat((nz-1:-1:0)', 1, nx)));
dmp = max(0, sp - id)/sp;
Vlist = [15 5 3];
for V = Vlist
  dzp = 2e-3; Tp = dzp/V; zs = 18e-3 + ((0:n-1) - (n-1)/2)*dzp;
  nsub = ceil(tRec/(0.6*dx/(cs*sqrt(2)))); dt = tRec/nsub; damp = exp(-dmp.^2*4*cs*dt/(sp*dx)*3);
  u = zeros(nz, nx); v = u; fr = zeros(nz, nx, nRec); t = 0; mu = rho*cs^2;
  G = zeros(nz, nx, n); for k = 1:n, G(:,:,k) = exp(-XX.^2/(2*(0.3e-3)^2) - (ZZ - zs(k)).^2/(2*(1.2e-3)^2)); end
  for r = 1:nRec
    for s = 1:nsub
      k = floor(t/Tp) + 1; f = zeros(nz, nx); if k <= n, f = G(:,:,k); end
      L = ([u(:,2:end) zeros(nz,1)] + [zeros(nz,1) u(:,1:end-1)] + [u(2:end,:); zeros(1,nx)] + [zeros(1,nx); u(1:end-1,:)] - 4*u)/dx^2;
      v = (v + dt/rho*(mu*L + 1e4*f)).*damp; u = u + dt*v; t = t + dt;
    end
    fr(:,:,r) = u;
  end
  M = V/cs; phi = asin(1/M); zrow = min(30e-3, 18e-3 + 5.5e-3/sqrt(M^2 - 1)); jr = round((zrow - z0)/dx) + 1;
  cols = (61 + 15):2:(61 + 40); tp = zeros(size(cols));
  for q = 1:numel(cols)
    vel = squeeze(diff(fr(jr, cols(q), :))); [~, im] = max(vel); d = 0;
    if im > 1 && im < numel(vel), a = vel(im-1); b = vel(im); g = vel(im+1); den = a - 2*b + g; if den < 0, d = 0.5*(a - g)/den; end, end
    tp(q) = (im + d - 0.5)*tRec;
  end
  p = polyfit(X(cols), tp, 1);
  fprintf('V = %4.1f m/s: M = %.2f, phi = %.1f deg, row z = %.1f mm, lateral TTP speed = %.2f m/s, c_s/cos(phi) = %.2f m/s\n', V, M, phi*180/pi, zrow*1e3, 1/p(1), cs/cos(phi));
  if V == Vlist(1), figure(1); imagesc(X*1e3, Zc*1e3, fr(:,:,31) - fr(:,:,30)); axis image; colormap(gray); xlabel('x (mm)'); ylabel('z (mm)'); title(sprintf('particle velocity at 3 ms, M = %.1f', M)); end
end
