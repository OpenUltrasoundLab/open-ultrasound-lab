%% m30_arfi_displacement.m
% Purpose : ARFI displacement from a radiation-force push. Axisymmetric (r,z)
%           finite-difference solution of the scalar Kelvin-Voigt shear-wave
%           equation rho*u_tt = mu*lap(u) + eta*lap(u_t) + F, Gaussian force
%           column (sigma_r = 0.35 mm, sigma_z = 2.5 mm). The step response is
%           computed once; a push of duration T is s(t) - s(t-T) (linearity).
%           F = 2*alpha*I/c sets the amplitude.
% Change  : cs (shear speed), eta, I (W/cm^2), a (dB/cm/MHz), f0, push list.
% Expect  : alpha = 23.0 Np/m, F = 1.50e5 N/m^3 at 500 W/cm^2; peak on-axis
%           displacement about 3.5 um for a 200 us push at 1.5 m/s (end-of-
%           push value F*T^2/(2 rho) = 3.0 um); long pushes approach the static
%           limit; time-to-peak at 2 and 4 mm gives back c_s; heating 12 mK.
rng(1);
rho = 1000; c = 1540;
I = 500; a = 0.5; f0 = 4;                       % W/cm^2, dB/cm/MHz, MHz
alpha = a*f0/8.686*100;                         % Np/m
F = 2*alpha*I*1e4/c;                            % N/m^3
fprintf('alpha = %.2f Np/m, F = 2 alpha I / c = %.3g N/m^3\n', alpha, F);
pushes = [50 100 200 400 800 1000]*1e-6;
dts = 1e-5; tEnd = 7e-3; ns = round(tEnd/dts) + 1;
dr = 1e-4; dz = 2e-4; nr = 81; nz = 101; sr = 0.35e-3; sz = 2.5e-3;
r = (0:nr-1)'*dr; z = ((0:nz-1) - (nz-1)/2)*dz;
f = exp(-r.^2/(2*sr^2)) * exp(-z.^2/(2*sz^2));   % nr x nz unit force
% absorbing layers
dd = zeros(nr, nz); spR = 20; spZ = 15;
for i = 1:nr, for j = 1:nz
  d = 0; if i > nr - spR, d = d + ((i - (nr - spR))/spR)^2; end
  jj = min(j - 1, nz - j); if jj < spZ, d = d + ((spZ - jj)/spZ)^2; end
  dd(i,j) = d;
end, end
lap = @(A) [4*(A(2,:) - A(1,:))/dr^2; ...
            (A(3:end,:) - 2*A(2:end-1,:) + A(1:end-2,:))/dr^2 + (A(3:end,:) - A(1:end-2,:))./(2*r(2:end-1)*dr); ...
            (-2*A(end,:) + A(end-1,:))/dr^2] + ...
           ([zeros(nr,1) A(:,1:end-1)] - 2*A + [A(:,2:end) zeros(nr,1)])/dz^2;
csList = [1.5 3]; eta = 0.5;
for ic = 1:numel(csList)
  cs = csList(ic); mu = rho*cs^2;
  lam = 8/dr^2 + 4/dz^2; dt = 0.7*2/(cs*sqrt(lam)); if eta > 0, dt = min(dt, 0.7*2*rho/(eta*lam)); end
  damp = exp(-dd*0.6*cs*dt/(spR*dr)*6);
  u = zeros(nr, nz); v = u; S = zeros(3, ns); jc = (nz+1)/2; t = 0; nxt = 1;
  while nxt <= ns
    while nxt <= ns && (nxt-1)*dts <= t + 1e-12, S(:,nxt) = [u(1,jc); u(21,jc); u(41,jc)]; nxt = nxt + 1; end
    v = (v + dt/rho*(mu*lap(u) + eta*lap(v) + f)).*damp; u = u + dt*v; t = t + dt;
  end
  % static limit (30.7): (1/(4 pi mu)) * int f/|r| dV
  [RR, ZZ] = meshgrid(((1:60) - 0.5)*0.05e-3, ((1:120) - 60.5)*0.1e-3);
  Is = sum(sum(exp(-RR.^2/(2*sr^2) - ZZ.^2/(2*sz^2))./sqrt(RR.^2 + ZZ.^2).*2*pi.*RR))*0.05e-3*0.1e-3;
  ust = F*Is/(4*pi*mu);
  fprintf('\nc_s = %.1f m/s, eta = %.1f Pa s (dt = %.2f us): static limit %.2f um\n', cs, eta, dt*1e6, ust*1e6);
  tt = (0:ns-1)*dts;
  for T = pushes
    k = round(T/dts); y = S(1,:) - [zeros(1,k) S(1,1:end-k)];
    [pk, ip] = max(y); irec = find(y(ip:end) < 0.5*pk, 1) + ip - 1;
    fprintf('  push %4.0f us: peak u = %6.2f um at %.2f ms, F T^2/(2 rho) = %6.2f um, recovery to 50%% %.2f ms after peak\n', ...
      T*1e6, F*pk*1e6, tt(ip)*1e3, F*T^2/(2*rho)*1e6, (irec - ip)*dts*1e3);
  end
  k = round(200e-6/dts); Y = S - [zeros(3,k) S(:,1:end-k)];
  tp = zeros(1,3); for q = 1:3, [~, ip] = max(Y(q,:)); de = Y(q,ip-1) - 2*Y(q,ip) + Y(q,ip+1); tp(q) = (ip - 1 + 0.5*(Y(q,ip-1) - Y(q,ip+1))/de)*dts; end
  fprintf('  TTP (200 us push) at r = 2 and 4 mm: %.2f and %.2f ms -> speed %.2f m/s\n', tp(2)*1e3, tp(3)*1e3, 2e-3/(tp(3) - tp(2)));
  if ic == 1, figure(1); plot(tt*1e3, F*Y'*1e6); xlabel('time (ms)'); ylabel('u (\mum)'); legend('r = 0', 'r = 2 mm', 'r = 4 mm'); title('200 \mus push, c_s = 1.5 m/s'); end
end
fprintf('\nHeating per 200 us push (no conduction): %.1f mK\n', 2*alpha*I*1e4*200e-6/3.7e6*1e3);
ap = a/8.686; zf = 4;                        % Np/cm/MHz, cm
fprintf('Optimum push frequency at %d cm: %.2f MHz (fixed source), %.2f MHz (focal gain f^2)\n', zf, 1/(2*ap*zf), 3/(2*ap*zf));
