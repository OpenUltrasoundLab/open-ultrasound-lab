%% m29_color_flow_wallfilter.m
% Purpose : Colour-flow processing chain on a simulated 2D ensemble:
%           clutter (40 dB above blood, slowly moving tissue) + blood in an
%           inclined vessel (parabolic profile) + white noise ->
%           polynomial-regression wall filter of degree p -> Kasai velocity,
%           power and normalized variance (1 - |R(1)|/R(0)) per pixel ->
%           power threshold. Also the filter's measured -3 dB cutoff,
%           compared with (0.44 + 0.41 p) PRF/N from Module 09.
% Change  : p (wall filter degree), N, PRF, CBR, tissue velocity vt.
% Expect  : N = 12, PRF = 2 kHz, f0 = 5 MHz: v_Nyq = 0.154 m/s. p = 0 leaves
%           tissue moving at 3 mm/s partly unfiltered (residual clutter >> blood,
%           velocities pulled toward the tissue value, flash outside the vessel);
%           p = 1-2 is the usable range; p = 3 eats slow flow near the walls.
rng(1);

%% Parameters
c = 1540; f0 = 5e6; PRF = 2e3; T = 1/PRF; N = 12;
theta = 60*pi/180;          % beam (vertical) to vessel angle
a = 3e-3;  vmax = 0.25;     % vessel radius (m), centreline speed (m/s)
CBRdB = 40;  SNRdB = 15;    % clutter-to-blood and blood-to-noise power ratios (dB)
vt = 3e-3;                  % axial tissue velocity (m/s), toward the probe
sigb = 40;  sigt = 3;       % blood and tissue spectral std (Hz)
nz = 60; nx = 100; dx = 0.4e-3;      % pixels (depth x lateral), pixel size (m)
vNyq = c*PRF/(4*f0);
fprintf('v_Nyq = %.3f m/s; axial centreline v = %.3f m/s; tissue Doppler = %.1f Hz\n', vNyq, vmax*cos(theta), 2*vt*f0/c);

[X, Z] = meshgrid(((1:nx) - nx/2)*dx, ((1:nz) - nz/2)*dx);
ang = pi/2 - theta;                          % vessel inclination from horizontal
dist = -X*sin(ang) + Z*cos(ang);             % signed distance from the vessel axis
lumen = abs(dist) < a;
vax = zeros(nz, nx);  vax(lumen) = vmax*(1 - (dist(lumen)/a).^2)*cos(theta);   % axial speed toward probe
fb  = 2*vax*f0/c;  ft = 2*vt*f0/c;
P = nz*nx;  n = (0:N-1);

% Signals (P pixels x N samples): sums of Kp phasors give Gaussian-like spectra
Kp = 8;
blood = zeros(P, N);  tiss = zeros(P, N);
for k = 1:Kp
    blood = blood + (randn(P,1) + 1j*randn(P,1))/sqrt(2*Kp) .* exp(2j*pi*(fb(:) + sigb*randn(P,1))*n*T);
    tiss  = tiss  + (randn(P,1) + 1j*randn(P,1))/sqrt(2*Kp) .* exp(2j*pi*(ft + sigt*randn(P,1))*n*T);
end
blood(~lumen(:), :) = 0;
noise = (randn(P,N) + 1j*randn(P,N))/sqrt(2) * 10^(-SNRdB/20);
Xsig = 10^(CBRdB/20)*tiss + blood + noise;

%% Regression filters, responses and processing
fax = linspace(-PRF/2, PRF/2, 801);  E = exp(2j*pi*n'*fax*T);        % N x F test tones
fprintf(' p  -3dB cutoff (Hz)  approx(0.44+0.41p)PRF/N  v_min(m/s)  resid CBR(dB)  detect%%  flash%%  mean|err| (m/s)\n');
maps = cell(1, 4);
for p = 0:3
    V = bsxfun(@power, ((0:N-1)' - (N-1)/2)/N, 0:p);  [Q, ~] = qr(V, 0);   % orthonormal polynomial basis
    A = eye(N) - Q*Q';                                   % regression (projection) filter
    H2 = sum(abs(A*E).^2, 1)/N;                          % |H(f)|^2 as in Module 09
    pos = fax >= 0;  fp = fax(pos);  hp = H2(pos);
    fc3 = fp(find(hp >= 0.5, 1));                        % first frequency reaching -3 dB
    Y = Xsig * A.';                                      % filter every pixel (rows)
    R0 = (sum(abs(Y(:,2:end)).^2, 2) + sum(abs(Y(:,1:end-1)).^2, 2))/2;
    R1 = sum(Y(:,2:end).*conj(Y(:,1:end-1)), 2);
    vk = c*PRF/(4*pi*f0)*angle(R1);
    nfl = 10^(-SNRdB/10) * (N - p - 1)/N * (N-1);       % expected noise-only R0
    dmask = R0 > 2*nfl;                                    % 3 dB above noise floor
    clutRes = mean(sum(abs(10^(CBRdB/20)*tiss*A.').^2, 2));  bl = mean(sum(abs(blood(lumen(:),:)*A.').^2, 2));
    inl = lumen(:) & dmask;
    fprintf(' %d     %6.0f             %6.0f               %.4f      %6.1f      %5.1f   %5.1f    %.4f\n', p, fc3, (0.44 + 0.41*p)*PRF/N, ...
        c*fc3/(2*f0*cos(theta)), 10*log10(clutRes/bl), 100*mean(dmask(lumen(:))), 100*mean(dmask(~lumen(:))), mean(abs(vk(inl) - vax(inl))));
    vm = vk;  vm(~dmask) = -1.1*vNyq;  maps{p+1} = reshape(vm, nz, nx);   % below colour range: black
    if p == 2, Hp2 = H2; end
end
fprintf('Detection threshold: R0 > 2 x expected noise power. Velocities in m/s are axial (toward probe > 0).\n');
fprintf('Slowest detectable flow (p = 2): about %.3f m/s along the vessel at theta = 60 deg\n', c*(0.44+0.82)*PRF/N/(2*f0*cos(theta)));

% Normalized variance map (p = 2): high where noise dominates
V = bsxfun(@power, ((0:N-1)' - (N-1)/2)/N, 0:2);  [Q, ~] = qr(V, 0);  A = eye(N) - Q*Q';  Y = Xsig*A.';
R0 = (sum(abs(Y(:,2:end)).^2, 2) + sum(abs(Y(:,1:end-1)).^2, 2))/2;  R1 = sum(Y(:,2:end).*conj(Y(:,1:end-1)), 2);
nv = 1 - abs(R1)./R0;
fprintf('Normalized variance 1-|R1|/R0 (p = 2): lumen centre %.2f, lumen edge %.2f, outside (noise) %.2f\n', ...
    median(nv(abs(dist(:)) < a/3)), median(nv(abs(dist(:)) > 0.75*a & lumen(:))), median(nv(~lumen(:))));

%% Figures
figure('Position', [40 40 1150 600]);
for p = 0:3
    subplot(2,3,p+1); imagesc(X(1,:)*1e3, Z(:,1)*1e3, maps{p+1}, [-1.1*vNyq vNyq]); axis image; colorbar;
    title(sprintf('Kasai v (m/s), regression p = %d', p)); xlabel('x (mm)'); ylabel('z (mm)');
end
colormap([0 0 0; jet(126)]);   % first entry (black) = below threshold
subplot(2,3,5); plot(fax, 10*log10(Hp2 + 1e-12), 'LineWidth', 1.5); hold on; plot(ft*[1 1], [-80 5], 'k--');
xlabel('Doppler frequency (Hz)'); ylabel('|H(f)|^2 (dB)'); title('Regression p = 2, N = 12'); ylim([-80 5]); grid on;
subplot(2,3,6); imagesc(X(1,:)*1e3, Z(:,1)*1e3, reshape(nv, nz, nx), [-0.01 1]); axis image; colorbar;
title('1 - |R(1)|/R(0), p = 2'); xlabel('x (mm)'); ylabel('z (mm)');
