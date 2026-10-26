%% m24_fk_migration_zero_offset.m  --  exploding reflector and Stolt migration of zero-offset data
% Purpose : zero-offset (monostatic synthetic aperture) data: every element fires alone and
%           records only its own echo, so the two-way delay to a point is 2r/c. That is the
%           one-way delay of a source at the point radiating at speed c/2 (the exploding
%           reflector). Stolt migration with v = c/2 maps the data spectrum along
%           f = (c/2) sqrt(Kx^2 + Kz^2) with Jacobian (c/2) Kz/sqrt(Kx^2 + Kz^2).
%           The script checks (1) Stolt against delay-and-sum (SAFT) on the same pixels,
%           (2) what using v = c instead of c/2 does, (3) the lambda/4 pitch rule for
%           zero-offset data: the data occupy |kx| <= 2f/c, so the pitch must satisfy
%           p <= c/(4 f_max) to avoid spatial aliasing at all angles.
% Change  : pitch list, v factor, target positions, f0.
% Expect  : Stolt and SAFT peaks at the true positions (errors of a few um), the same
%           lateral width (two-way, full 38.4 mm aperture: ~0.1-0.2 mm) within a few %.
%           v = c puts the point at twice its depth and defocuses it. Artifact level away
%           from the targets drops sharply from pitch 0.30 mm (lambda) to 0.075 mm (lambda/4).
% Model   : linear acoustics, Born, point-like elements, no attenuation or spreading,
%           Gaussian pulse, IQ at 20 MHz. Runtime ~10 s in Octave.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);
tStart = tic;

%% Parameters
c = 1.540; f0 = 5; lam = c/f0; FBW = 0.6; fs = 20;   % mm/us, MHz, mm, -, MHz
L = 38.4;                                            % aperture (mm), same for every pitch
pitches = [0.30 0.15 0.075];                         % lambda, lambda/2, lambda/4 at 5 MHz
ptsX = [0 -6 6 0].';  ptsZ = [10 20 20 30].';
nt = 1100;  ntFFT = 4096;  tc = nt/(2*fs);           % samples, FFT length, centering (us)
sig_f = FBW*f0/(2*sqrt(2*log(2)));  sig_t = 1/(2*pi*sig_f);
kt = (-ceil(4*sig_t*fs):ceil(4*sig_t*fs)).';  envK = exp(-(kt/fs).^2/(2*sig_t^2));
dxf = 0.075;                                         % common lateral image grid (mm)
Fnum = 1.5;  alphaMax = atan(1/(2*Fnum));           % receive F# (SAFT aperture, Stolt fan)
NZ = 1024;  dKz = 2*(fs/2048)/c;  Kz = (0:NZ-1).'*dKz;  zf = (0:NZ-1).'/(NZ*dKz);
Dfq = fs/ntFFT;  fq = [0:ntFFT/2-1, -ntFFT/2:-1]*Dfq + f0;
bq = find(fq >= 0.8 & fq <= 9.2);  [fB, o] = sort(fq(bq));  bq = bq(o);  fB = fB(:);  nB = numel(fB);
fprintf('zero-offset aliasing-free pitch c/(4 f) = %.3f mm at f0, %.3f mm at %.1f MHz (-24 dB band edge)\n\n', ...
    c/(4*f0), c/(4*8), 8);

fprintf('%-6s %-6s %-6s %8s %8s %9s %9s %10s\n', 'pitch', 'p/lam', 'method', 'target', 'lat-6dB', 'dx peak', 'dz peak', 'artifact');
fprintf('%-6s %-6s %-6s %8s %8s %9s %9s %10s\n', '(mm)', '', '', '(x,z)', '(mm)', '(um)', '(um)', '(dB)');
store = cell(1, numel(pitches));
for ip = 1:numel(pitches)
    p = pitches(ip);  Nel = round(L/p);  xe = ((0:Nel-1) - (Nel-1)/2)*p;
    %% zero-offset channel data, IQ (nt x Nel)
    tau = 2*sqrt((ptsX - xe).^2 + ptsZ.^2)/c;            % two-way delay (us)
    w = exp(-1i*2*pi*f0*tau);  q = tau*fs;  i0 = floor(q);  fr = q - i0;
    li = i0 + 1 + nt*(0:Nel-1);
    Sd = accumarray(li(:), w(:).*(1-fr(:)), [nt*Nel 1]) + accumarray(li(:)+1, w(:).*fr(:), [nt*Nel 1]);
    S  = conv2(reshape(Sd, nt, Nel), envK, 'same');
    %% Stolt (v = c/2), lateral zero-padding to 4x the aperture, upsampled to dxf
    nxp = 2^nextpow2(4*Nel);  U = round(p/dxf);  nxf = U*nxp;
    xf  = xe(1) + (0:nxf-1)*dxf;  Kx = [0:nxf/2-1, -nxf/2:-1]/(nxf*dxf);
    [KX, KZ] = meshgrid(Kx, Kz);
    D  = fft(S, ntFFT, 1);  D = D(bq, :) .* exp(1i*2*pi*fB*tc);  D = fft(D, nxp, 2);
    kxd = [0:nxp/2-1, -nxp/2:-1]/(nxp*p);  Dpre = D;
    D  = D .* (abs(kxd) <= 2*fB*sin(alphaMax)/c);       % fan filter: zero-offset kx = 2 f sin(alpha)/c
    Du = zeros(nB, nxf);  Du(:, 1:nxp/2) = D(:, 1:nxp/2);  Du(:, end-nxp/2+1:end) = D(:, nxp/2+1:end);
    imgs = cell(1,2);
    for vf = [0.5 1]                                     % velocity c/2 (correct) and c (wrong)
        v = vf*c;  F = v*sqrt(KX.^2 + KZ.^2);  J = v*KZ./(sqrt(KX.^2 + KZ.^2) + eps);
        qq = (F - fB(1))/Dfq + 1;  k0 = floor(qq);  ww = qq - k0;  ok = k0 >= 1 & k0 < nB;
        lk = k0(ok) + nB*(ceil(find(ok)/NZ) - 1);
        I = zeros(NZ, nxf);  I(ok) = ((1 - ww(ok)).*Du(lk) + ww(ok).*Du(lk+1)) .* J(ok) .* exp(-1i*2*pi*F(ok)*tc);
        imgs{1 + (vf == 1)} = ifft2(I);
    end
    E = abs(imgs{1});  E = E/max(E(:));
    store{ip} = struct('E', E, 'xf', xf, 'Ew', abs(imgs{2})/max(abs(imgs{1}(:))), 'S', S, 'xe', xe, 'D', Dpre);
    %% SAFT (delay and sum, dynamic aperture F#) on the Stolt pixels around each target
    for it = 1:numel(ptsX)
        jx = find(abs(xf - ptsX(it)) <= 3);  jz = find(abs(zf - ptsZ(it)) <= 1);
        [XL, ZL] = meshgrid(xf(jx), zf(jz));
        Td = 2*sqrt((XL(:) - xe).^2 + ZL(:).^2)/c;       % pixels x elements (us)
        qd = Td*fs;  kd = floor(qd);  fd = qd - kd;
        okd = kd >= 0 & kd < nt - 1 & abs(XL(:) - xe) <= ZL(:)/(2*Fnum);
        kd(~okd) = 0;  ld = kd + 1 + nt*(0:Nel-1);
        V  = ((1-fd).*S(ld) + fd.*S(ld+1)) .* okd .* exp(1i*2*pi*f0*Td);
        Bd = reshape(abs(sum(V, 2)), size(XL));
        Bs = E(jz, jx);
        for m = 1:2
            if m == 1, B = Bs; nm = 'Stolt'; else, B = Bd; nm = 'SAFT'; end
            B = B/max(B(:));  [~, k] = max(B(:));  [kz, kx] = ind2sub(size(B), k);
            xs = XL(1,:);  zs = ZL(:,1);  xq = xs(kx) + (-0.6:0.002:0.6);  zq = zs(kz) + (-0.4:0.002:0.4);
            Lc = interp1(zs, B(:, kx), zq, 'spline');  [~, b2] = max(Lc);  zpk = zq(b2);
            rowv = interp1(zs, B, zpk, 'spline');  Lr = interp1(xs, rowv, xq, 'spline');
            [pr, a2] = max(Lr);  xpk = xq(a2);  wl = sum(Lr >= pr/2)*0.002;
            far = abs(XL - xpk) > 1.0 | abs(ZL - zpk) > 0.5;
            fprintf('%-6.3f %-6.2f %-6s %8s %8.3f %9.1f %9.1f %10.1f\n', p, p/lam, nm, sprintf('%g,%g', ptsX(it), ptsZ(it)), ...
                wl, (xpk - ptsX(it))*1e3, (zpk - ptsZ(it))*1e3, 20*log10(max(B(far))));
        end
    end
    %% wrong velocity: where does the target at 10 mm go?
    Ew = store{ip}.Ew;  sel = abs(xf) < 3;  [pw, k] = max(max(Ew(zf < 26, sel), [], 2));
    fprintf('       v = c instead of c/2: brightest point near x = 0 at z = %.2f mm (true 10), peak %.1f dB re correct\n\n', zf(k), 20*log10(pw));
end
fprintf('uniform kx support |kx| <= 2 sin(alpha)/lam gives a -6 dB width 1.207 lam/(4 sin alpha) = %.3f mm\n', 1.207*lam/(4*sin(alphaMax)));
fprintf('fan edge reaches the pitch Nyquist 1/(2p) at f = c/(4 p sin alpha): %.1f MHz (p = 0.30), %.1f MHz (0.15), %.1f MHz (0.075)\n', ...
    c./(4*pitches*sin(alphaMax)));
fprintf('total run time %.1f s\n', toc(tStart));

%% Figures
s = store{1};  t = (0:nt-1)/fs;
figure('Position', [50 50 1150 650]);
subplot(2,3,1); imagesc(s.xe, t, real(s.S)); colormap(gray); xlabel('element x (mm)'); ylabel('t (us)');
title('zero-offset data (real IQ), pitch 0.30 mm');
subplot(2,3,2); nxp = size(s.D,2); kxd = [0:nxp/2-1, -nxp/2:-1]/(nxp*pitches(1));
P = 20*log10(abs(fftshift(s.D, 2)) + eps);  P = P - max(P(:));
imagesc(fftshift(kxd), fB, max(P, -50)); axis xy; hold on; axis([-1.8 1.8 fB(1) fB(end)]);
plot(2*fB/c, fB, 'r--', -2*fB/c, fB, 'r--'); plot([1 1]/(2*pitches(1)), fB([1 end]), 'y', -[1 1]/(2*pitches(1)), fB([1 end]), 'y');
xlabel('k_x (cycles/mm)'); ylabel('f (MHz)'); title('|D(k_x,f)| (dB), p = 0.30 mm: red |k_x| = 2f/c (folds), yellow Nyquist');
subplot(2,3,3); s = store{1}; imagesc(s.xf, zf, max(20*log10(s.E + eps), -50), [-50 0]); axis image; axis([-10 10 5 35]);
title('Stolt, v = c/2, pitch 0.30 mm (dB)'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,4); s = store{3}; imagesc(s.xf, zf, max(20*log10(s.E + eps), -50), [-50 0]); axis image; axis([-10 10 5 35]);
title('Stolt, v = c/2, pitch 0.075 mm (dB)'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,5); s = store{3}; imagesc(s.xf, zf, max(20*log10(s.Ew/max(s.Ew(:)) + eps), -50), [-50 0]); axis image; axis([-15 15 5 70]);
title('Stolt with v = c (wrong): depth x2, defocused'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(2,3,6); s1 = store{1}; s3 = store{3}; [~, jz] = min(abs(zf - 20));
plot(s1.xf, 20*log10(s1.E(jz,:) + eps), 'r', s3.xf, 20*log10(s3.E(jz,:) + eps), 'k'); axis([-12 12 -60 0]); grid on;
xlabel('x (mm)'); ylabel('dB'); legend('pitch 0.30 mm', 'pitch 0.075 mm'); title('lateral profile at z = 20 mm');
