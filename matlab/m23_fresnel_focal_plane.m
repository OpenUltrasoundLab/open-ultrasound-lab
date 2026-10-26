%% m23_fresnel_focal_plane.m  --  Fresnel vs Fraunhofer: why the focus is a Fourier plane
% Purpose : CW field of a 1D aperture (a line of point sources, 1/r spreading, 3D Green's
%           function), computed by direct summation without the paraxial approximation.
%           (1) Unfocused: on-axis intensity and its last maximum vs D^2/(4 lambda).
%           (2) Focused at F: lateral profile at z = F vs the Fourier transform of the
%               apodization (rect -> sinc, Hann), -6 dB width and first sidelobe.
%           (3) Focal shift and depth of field vs Fresnel number N_F = D^2/(4 lambda F):
%               the intensity peak moves toward the array and the -3 dB depth of field
%               departs from 6.95 lambda F#^2 when N_F is small.
% Change  : f0, D, F, apodization ('rect' or 'hann').
% Expect  : last axial maximum at ~1.0 D^2/(4 lambda) = 299 mm (D = 19.2 mm, 5 MHz).
%           At F = 30 mm (F# 1.56, N_F = 10) the direct field is 3 % wider than the paraxial
%           Fourier prediction and its first sidelobe is -12.7 instead of -13.3 dB (rect):
%           the residual is the non-paraxial error, and it shrinks at F = 60 mm (F# 3.1).
%           Peak/F = 0.99 at N_F = 10, 0.96 at N_F = 5, ~0.90 at N_F = 3, 0.6 at N_F = 1.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

%% Parameters
c   = 1.540;            % mm/us
f0  = 5;                % MHz
lam = c/f0;             % mm
k   = 2*pi/lam;         % rad/mm
D   = 19.2;             % aperture length (mm) = 64 elements x 0.3 mm
F   = 30;               % focal depth (mm)
apod = 'rect';          % 'rect' or 'hann'
nSrc = 400;             % point sources across the aperture (spacing 0.048 mm < lambda/6)
xi  = linspace(-D/2, D/2, nSrc);                     % source positions (mm), row
u   = 2*xi/D;
if strcmp(apod, 'hann'), a = 0.5 + 0.5*cos(pi*u); else, a = ones(size(u)); end
field = @(x, z, foc) sum( (ones(numel(x),1)*a) .* exp(-1i*k*(sqrt((x(:)-xi).^2 + z^2) - foc)) ...
                     ./ sqrt((x(:)-xi).^2 + z^2), 2);   % complex CW field on a lateral line
focLaw = @(Fz) sqrt(xi.^2 + Fz^2) - Fz;              % focusing path correction (mm)

%% (1) Unfocused on-axis intensity
zAx = linspace(2, 600, 3000);
I0  = zeros(size(zAx));
for iz = 1:numel(zAx), I0(iz) = abs(field(0, zAx(iz), 0*xi))^2; end
isMax = [false, I0(2:end-1) > I0(1:end-2) & I0(2:end-1) > I0(3:end), false];
zLast = max(zAx(isMax));
fprintf('(1) unfocused: last axial maximum %.0f mm; D^2/(4 lam) = %.0f mm (ratio %.2f)\n', ...
    zLast, D^2/(4*lam), zLast/(D^2/(4*lam)));
fprintf('    Fraunhofer distances: D^2/lam = %.0f mm (edge phase pi/4), 2D^2/lam = %.0f mm (pi/8)\n', ...
    D^2/lam, 2*D^2/lam);

%% (2) Focal-plane lateral profile vs Fourier transform of the apodization
x   = linspace(-8, 8, 8001);  mid = (numel(x)+1)/2;
w6  = @(p) 2*x(find(p >= 0.5, 1, 'last'));            % symmetric profile, full width
for Fz = [F 2*F]
    p = abs(field(x, Fz, focLaw(Fz)));  p = p/max(p);
    % Fourier prediction: A(x/(lam Fz)) for the same apodization (paraxial)
    Af = abs(exp(-1i*2*pi*(x(:)/(lam*Fz))*xi) * a.');  Af = Af/max(Af);
    iN = mid - 1 + find(diff(p(mid:end)) > 0, 1);      % first minimum right of center
    iA = mid - 1 + find(diff(Af(mid:end)) > 0, 1);
    fprintf('(2) focused at F = %g mm (F# = %.2f, N_F = %.1f), %s apodization:\n', Fz, Fz/D, D^2/(4*lam*Fz), apod);
    fprintf('    -6 dB width: field %.3f mm, Fourier %.3f mm (ratio %.3f); 1.207 lam F# = %.3f mm\n', ...
        w6(p), w6(Af), w6(p)/w6(Af), 1.207*lam*Fz/D);
    fprintf('    first sidelobe: field %.2f dB, Fourier %.2f dB\n', 20*log10(max(p(iN:end))), 20*log10(max(Af(iA:end))));
    if Fz == F, pF = p; A = Af; xF = x; end
end
x = xF;

%% (3) Focal shift and depth of field vs Fresnel number
Fs = [20 30 40 60 80 100 150 200 300];
fprintf('(3) %6s %6s %6s %8s %10s %12s\n', 'F(mm)', 'F#', 'N_F', 'peak/F', 'DOF(mm)', '6.95lamF#^2');
tab = zeros(numel(Fs), 4);
for j = 1:numel(Fs)
    Fz = Fs(j);  zz = linspace(0.35*Fz, 2.5*Fz, 2500);   Iz = zeros(size(zz));   % start past the 1/r near-field rise
    fl = focLaw(Fz);
    for iz = 1:numel(zz), Iz(iz) = abs(field(0, zz(iz), fl))^2; end
    [Im, ip] = max(Iz); lo = ip; hi = ip;
    while lo > 1 && Iz(lo) >= Im/2, lo = lo - 1; end
    while hi < numel(zz) && Iz(hi) >= Im/2, hi = hi + 1; end
    NF = D^2/(4*lam*Fz);  dof = zz(hi) - zz(lo);
    tab(j,:) = [NF, zz(ip)/Fz, dof, 6.95*lam*(Fz/D)^2];
    fprintf('    %6g %6.2f %6.2f %8.3f %10.1f %12.1f\n', Fz, Fz/D, NF, zz(ip)/Fz, dof, 6.95*lam*(Fz/D)^2);
end

%% Field maps (unfocused and focused at F) and profiles
xm = linspace(-8, 8, 161); zm = linspace(2, 80, 196);
M0 = zeros(numel(zm), numel(xm)); M1 = M0;
for iz = 1:numel(zm)
    M0(iz,:) = abs(field(xm, zm(iz), 0*xi)).';
    M1(iz,:) = abs(field(xm, zm(iz), focLaw(F))).';
end
figure('Position', [40 40 1150 760]);
subplot(2,3,1); imagesc(xm, zm, 20*log10(M0/max(M0(:))), [-40 0]); colormap(gray); colorbar;
title('unfocused CW field (dB)'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,2); imagesc(xm, zm, 20*log10(M1/max(M1(:))), [-40 0]); colorbar;
title(sprintf('focused at %g mm (dB)', F)); xlabel('x (mm)'); ylabel('z (mm)');
IdB = 10*log10(I0/interp1(zAx, I0, zLast));            % dB re the last axial maximum
subplot(2,3,3); plot(zAx, IdB, 'LineWidth', 1.2); hold on; xlim([20 600]); ylim([-12 8]);
plot([1 1]*D^2/(4*lam), [-12 8], 'k--'); grid on; xlabel('z (mm)'); ylabel('on-axis intensity (dB re last max)');
title('unfocused: last maximum near D^2/(4\lambda)');
subplot(2,3,4); plot(x, 20*log10(pF + 1e-6), x, 20*log10(A + 1e-6), '--', 'LineWidth', 1.2); grid on;
xlim([-4 4]); ylim([-50 3]); xlabel('x (mm) at z = F'); ylabel('dB'); legend('direct field', 'FT of apodization');
title('focal plane = Fourier plane');
subplot(2,3,5); semilogx(tab(:,1), tab(:,2), 'o-', 'LineWidth', 1.2); grid on;
xlabel('Fresnel number N_F = D^2/(4\lambda F)'); ylabel('intensity peak / F'); title('focal shift');
subplot(2,3,6); semilogx(tab(:,1), tab(:,3)./tab(:,4), 's-', 'LineWidth', 1.2); grid on;
xlabel('Fresnel number N_F'); ylabel('DOF / (6.95 \lambda F#^2)'); title('depth-of-field rule vs N_F');
