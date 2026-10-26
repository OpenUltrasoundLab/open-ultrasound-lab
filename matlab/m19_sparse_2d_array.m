%% m19_sparse_2d_array.m  --  2D matrix array: full, periodic and random sparse layouts
% Purpose : narrowband far-field pattern of a 32 x 32 lambda/2 matrix by 2D FFT, in sine space
%           (u, v) = (sin th cos phi, sin th sin phi). Compares the fully populated aperture, a
%           periodic layout (every 2nd element in x and y) and random layouts with K elements.
%           For random layouts the mean sidelobe power relative to the peak is predicted by the
%           random-walk formula (N^2 - K)/(K (N^2 - 1)) ~ 1/K.
% Change  : K values, number of random draws, pad size P; spiral taper wfun, aperture R.
% Expect  : full PSL -13.3 dB (rect, along the axes; -13.4 on this FFT grid); periodic every-2nd: grating lobes
%           at 0 dB at endfire (|u| = 1), and at u = -0.5 when steered to u0 = 0.5;
%           random K = 128 / 256 / 512 (20 draws): mean sidelobe about -21.6 / -25.0 / -29.2 dB vs
%           -21.6 / -25.3 / -30.1 dB predicted; peak sidelobe of a single draw about -9 to -15 dB.
%           Channel count of the full aperture: 1024; data rate 1024 x 40 MHz x 12 bit = 492 Gb/s.
%           Off-grid layouts, 256 elements on a 10 mm aperture at 5 MHz (direct sum): random on a 64 x 64
%           lambda/2 grid, Fermat spiral, density-tapered spiral: mean sidelobe about -24.4 / -24.8 / -24.4 dB
%           (1/K = -24.1 dB), PSL -15.3 / -13.8 / -15.1 dB, -6 dB main-lobe width 0.037 / 0.043 / 0.051 in u.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

N = 32; P = 256;                                     % elements per side, FFT size
fprintf('Full matrix: %d elements; raw data %.0f Gb/s at 40 MHz, 12 bit\n', N^2, N^2*40e6*12/1e9);
uu = ((0:P-1) - P/2) / P * 2;                        % u for d = lambda/2: u = (lambda/d) f = 2 f
[U, V] = meshgrid(uu, uu);
vis  = (U.^2 + V.^2) <= 1;                           % visible region
ml   = 2/N;                                          % first null of the full aperture (u = lambda/(N d))
main = (abs(U) < ml) & (abs(V) < ml);
side = vis & ~main;
pat  = @(mask) fftshift(abs(fft2(mask, P, P)));      % |array factor|, unnormalized

%% Full and periodic layouts
A = pat(ones(N)); A = A / max(A(:));
fprintf('full aperture        : PSL %6.1f dB, mean sidelobe %6.1f dB\n', 20*log10(max(A(side))), 10*log10(mean(A(side).^2)));
mp = zeros(N); mp(1:2:end, 1:2:end) = 1;
Ap = pat(mp); Ap = Ap / max(Ap(:));
fprintf('periodic every 2nd   : %d elements, highest lobe outside main %5.1f dB (grating lobes at endfire, |u| = 1)\n', sum(mp(:)), 20*log10(max(Ap(side))));
u0 = 0.5; ramp = repmat(exp(1i*pi*u0*(0:N-1)), N, 1);   % steer to u0 along x (d = lambda/2: phase pi u0 per element)
Aps = pat(mp .* ramp); Aps = Aps / max(Aps(:));
sides = vis & ~((abs(U - u0) < ml) & (abs(V) < ml));
[~, ig] = max(Aps(:) .* sides(:));
fprintf('periodic, steered u0 = %.1f: highest lobe %5.1f dB at u = %.2f (grating lobe at u0 - 1)\n', u0, 20*log10(Aps(ig)), U(ig));

%% Random layouts
Ks = [128 256 512]; nd = 20;
fprintf('\n   K   mean SL (sim)   predicted   PSL mean (min..max)\n');
for K = Ks
    msl = zeros(1, nd); psl = zeros(1, nd);
    for t = 1:nd
        m = zeros(N); idx = randperm(N^2); m(idx(1:K)) = 1;
        Ar = pat(m); Ar = Ar / max(Ar(:));
        msl(t) = 10*log10(mean(Ar(side).^2)); psl(t) = 20*log10(max(Ar(side)));
        if K == 256 && t == 1, m256 = m; A256 = Ar; end
    end
    fprintf('%5d     %6.1f        %6.1f      %5.1f (%5.1f..%5.1f)\n', K, mean(msl), 10*log10((N^2 - K)/(K*(N^2 - 1))), mean(psl), min(psl), max(psl));
end

%% Off-grid layouts with 256 elements on a 10 mm aperture (same layouts as the page drawing)
% Direct phasor sum (positions are not on a grid, so no FFT): P(u,v) = sum exp(j 2 pi (x u + y v)/lambda).
lam = 0.308; R = 5; Ke = 256; n = (0:Ke-1)'; ga = pi*(3 - sqrt(5));   % mm, golden angle 137.5 deg
rF = R*sqrt((n + 0.5)/Ke);                                           % Fermat spiral: uniform density
wfun = @(r) 0.25 + 0.75*0.5*(1 + cos(pi*r/R));                       % density taper (pedestal + raised cosine)
rr = linspace(0, R, 4001); cw = [0, cumsum(wfun((rr(1:end-1) + rr(2:end))/2) .* 2*pi.*(rr(1:end-1) + rr(2:end))/2 .* diff(rr))];
rT = interp1(cw, rr, (n + 0.5)/Ke*cw(end));                          % invert the cumulative count
idx = randperm(64^2, Ke)'; gx = (mod(idx - 1, 64) - 31.5)*0.154; gy = (floor((idx - 1)/64) - 31.5)*0.154;
lay = {gx, gy; rF.*cos(n*ga), rF.*sin(n*ga); rT.*cos(n*ga), rT.*sin(n*ga)};
names = {'random 256 of 64 x 64 (lambda/2)', 'Fermat spiral, 10 mm', 'density-tapered spiral'};
ug = linspace(-1, 1, 161); [Ug, Vg] = meshgrid(ug, ug); visg = (Ug.^2 + Vg.^2) <= 1;
uc = linspace(-0.15, 0.15, 601);                                      % fine cut for the main lobe
fprintf('\n10 mm aperture, 256 elements, 5 MHz (1/K = %.1f dB):\n', 10*log10(1/Ke));
fprintf('   layout                              PSL     mean SL   -6 dB width (u)\n');
Aoff = cell(1, 3);
for k = 1:3
    x = lay{k, 1}; y = lay{k, 2};
    A2 = abs(exp(1i*2*pi/lam*(Ug(:)*x' + Vg(:)*y')) * ones(Ke, 1)) / Ke; A2 = reshape(A2, size(Ug));
    ac = abs(exp(1i*2*pi/lam*(uc(:)*x')) * ones(Ke, 1)) / Ke;
    w6 = 2*max(uc(ac >= 0.5));                                        % full -6 dB width of the main lobe
    sideg = visg & (sqrt(Ug.^2 + Vg.^2) > 1.5*w6);                     % exclude the main lobe
    fprintf('   %-34s %6.1f    %6.1f      %.4f\n', names{k}, 20*log10(max(A2(sideg))), 10*log10(mean(A2(sideg).^2)), w6);
    Aoff{k} = A2;
end

%% Figures
dBimg = @(X) max(20*log10(X + 1e-12), -40) .* vis - 45*(~vis);
figure('Name', 'Sparse 2D arrays');
subplot(2, 3, 1); imagesc(ones(N)); axis image; title('full 32 x 32'); colormap(gray);
subplot(2, 3, 2); imagesc(mp); axis image; title('periodic every 2nd');
subplot(2, 3, 3); imagesc(m256); axis image; title('random K = 256');
subplot(2, 3, 4); imagesc(uu, uu, dBimg(A));    axis image; caxis([-40 0]); xlabel('u'); ylabel('v'); title('full (dB)');
subplot(2, 3, 5); imagesc(uu, uu, dBimg(Aps));  axis image; caxis([-40 0]); xlabel('u'); title('periodic, steered u_0 = 0.5');
subplot(2, 3, 6); imagesc(uu, uu, dBimg(A256)); axis image; caxis([-40 0]); xlabel('u'); title('random 256 (dB)'); colorbar;
figure('Name', 'Cuts v = 0');
plot(uu, 20*log10(A(P/2+1, :) + 1e-12), 'k', uu, 20*log10(A256(P/2+1, :) + 1e-12), 'b'); hold on;
plot([-1 1], 10*log10((N^2 - 256)/(256*(N^2 - 1)))*[1 1], 'r--');
ylim([-60 0]); grid on; xlabel('u = sin \theta (v = 0)'); ylabel('dB'); legend('full', 'random K = 256', '1/K floor');
figure('Name', 'Off-grid 256-element layouts on 10 mm');
dBg = @(X) max(20*log10(X + 1e-12), -40) .* visg - 45*(~visg);
for k = 1:3
    subplot(2, 3, k); plot(lay{k, 1}, lay{k, 2}, '.', 'MarkerSize', 6); axis image; axis([-5.2 5.2 -5.2 5.2]); xlabel('x (mm)'); ylabel('y (mm)'); title(names{k});
    subplot(2, 3, 3 + k); imagesc(ug, ug, dBg(Aoff{k})); axis image; caxis([-40 0]); xlabel('u'); ylabel('v'); title('|P(u,v)| (dB)');
end
colormap(gray); colorbar;
