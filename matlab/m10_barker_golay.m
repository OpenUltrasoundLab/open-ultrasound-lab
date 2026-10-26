%% m10_barker_golay.m
% Binary phase codes: Barker-13 (matched and least-squares mismatched filter)
% and Golay complementary pairs, including what motion between the two
% firings does to the Golay sidelobe cancellation.
%   - Barker-13: autocorrelation peak 13, sidelobes 0 or 1 -> PSL = -22.28 dB.
%   - LS mismatched filter of length L: minimizes the summed squared sidelobes;
%     PSL falls with L for about 0.2 dB of SNR loss.
%   - Golay pair (A,B): R_A + R_B = 2N delta. If echo B is phase-rotated by phi
%     (axial displacement d between firings, phi = 4 pi f0 d / c), the residual
%     PSL is (max|R_A(k)|/N) tan(phi/2).
%   - RF check: 1-cycle chips at f0, Gaussian transducer, true delay 2d/c.
% What to change: Lmm, Nlist, f0, fbw, d_um.
% Expected output: Barker PSL -22.28 dB, gain 11.14 dB; LS filter L = 39: PSL
% -38.5 dB, loss 0.20 dB; Golay N = 16: single-code PSL -10.1 dB, pair sum
% sidelobes exactly 0, d = 2 um -> -37.9 dB (phase-only theory), -38.4 dB on RF.
% Runtime about 1 s.
clear; close all; rng(1);

dB = @(x) 20*log10(max(abs(x), 1e-300));
c  = 1540;                 % sound speed (m/s)
f0 = 5e6;                  % centre frequency (Hz)

%% 1. Barker-13, matched filter
b13 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]';
r = conv(b13, flipud(b13));                      % autocorrelation, peak at index 13
side = r; side(13) = [];
fprintf('Barker-13 autocorrelation: %s\n', mat2str(r'));
fprintf('PSL = %.2f dB, compression gain 10log10(13) = %.2f dB\n', dB(max(abs(side))/13), 10*log10(13));

%% 2. Least-squares mismatched filters for Barker-13
Lmm = [13 21 27 39 51 65];
fprintf('\nLS mismatched filter (target: unit peak, zero sidelobes)\n    L    PSL (dB)   SNR loss (dB)\n');
hbest = [];
for L = Lmm
    M = 13 + L - 1; C = zeros(M, L);
    for j = 1:L, C(j:j+12, j) = b13; end           % convolution matrix: y = C*h
    d = zeros(M,1); d((M+1)/2) = 1;                % desired output: a centred spike
    h = C \ d; y = C*h;
    [pk, k] = max(abs(y)); ys = y; ys(k) = [];
    loss = 10*log10(pk^2/(sum(h.^2)*sum(b13.^2)));  % relative to the matched filter (13)
    fprintf('%5d   %7.1f     %6.2f\n', L, dB(max(abs(ys))/pk), loss);
    if L == 39, hbest = h; ybest = y/pk; end
end

%% 3. Golay complementary pairs
Nlist = [4 8 16 32 64];
fprintf('\nGolay pairs\n    N   max|R_A|   PSL single (dB)   max|R_A+R_B| off-peak   gain 10log10(2N)\n');
for N = Nlist
    a = 1; b = 1;
    while numel(a) < N, t = [a b]; b = [a -b]; a = t; end   % recursive construction
    ra = conv(a, fliplr(a)); rb = conv(b, fliplr(b)); s = ra + rb;
    sa = ra; sa(N) = []; ss = s; ss(N) = [];
    fprintf('%5d   %5d       %7.1f             %g                 %.2f dB\n', N, max(abs(sa)), dB(max(abs(sa))/N), max(abs(ss)), 10*log10(2*N));
end
N = 16; a = 1; b = 1;
while numel(a) < N, t = [a b]; b = [a -b]; a = t; end
ra = conv(a, fliplr(a)); rb = conv(b, fliplr(b)); sa = ra; sa(N) = []; mA = max(abs(sa));
d_um = [0.5 1 2 5 10 20 38.5];
fprintf('\nGolay N = 16, phase rotation phi = 4 pi f0 d / c between firings (f0 = %.0f MHz)\n   d (um)   phi (rad)   PSL simulated   PSL theory (mA/N) tan(phi/2)\n', f0/1e6);
for d = d_um
    phi = 4*pi*f0*d*1e-6/c; z = ra + exp(1j*phi)*rb; zs = z; zs(N) = [];
    fprintf('%8.1f   %7.3f      %7.1f dB        %7.1f dB\n', d, phi, dB(max(abs(zs))/abs(z(N))), dB(mA/N*tan(phi/2)));
end

%% 4. RF simulation: 1-cycle chips, transducer band limit, true delay 2d/c on firing B
fs = 16*f0; nc = round(fs/f0);                    % 16 samples per chip (one carrier cycle)
fbw = 0.7; sf = fbw*f0/2.3548;                    % two-way Gaussian transducer response
chip = sin(2*pi*(0:nc-1)'/nc);
codewave = @(code) reshape(chip*code(:)', [], 1);
Nfft = 4096; fax = [0:Nfft/2, -Nfft/2+1:-1]'*fs/Nfft;
Htr = exp(-(abs(fax)-f0).^2/(2*sf^2));
wA = real(ifft(fft(codewave(a), Nfft).*Htr)); wB = real(ifft(fft(codewave(b), Nfft).*Htr));   % received pulse shapes
an = zeros(Nfft,1); an(1) = 1; an(2:Nfft/2) = 2; an(Nfft/2+1) = 1;     % analytic-signal weights
env = @(Y) abs(ifft(Y.*an));
WA = fft(wA); WB = fft(wB);
fprintf('\nRF simulation (N = 16, 1 cycle/chip, %.0f%% transducer)\n   d (um)   PSL sum (dB)   phase-only theory (dB)\n', fbw*100);
dplot = [0 5 20]; envs = cell(1,3);
for d = [0 0.5 2 5 20]
    tau = 2*d*1e-6/c;                             % extra round-trip delay of echo B (s)
    YB = WB.*exp(-1j*2*pi*fax*tau);               % exact fractional delay (band-limited signal)
    Ysum = WA.*conj(WA) + YB.*conj(WB);           % each firing compressed by its own matched filter
    e = fftshift(env(Ysum)); [pk,k] = max(e); en = e/pk;
    % sidelobe region: more than 4 chips from the peak (the band-limited mainlobe skirt is excluded)
    mask = abs((1:Nfft)' - k) > 4*nc;
    psl = dB(max(en(mask)));
    phi = 4*pi*f0*d*1e-6/c;
    if d == 0, th = -Inf; else, th = dB(mA/N*tan(phi/2)); end
    fprintf('%8.1f      %7.1f          %7.1f\n', d, psl, th);
    ii = find(dplot == d); if ~isempty(ii), envs{ii} = en; end
end
eA = fftshift(env(WA.*conj(WA))); eA = eA/max(eA);

%% Figures
figure('Position',[100 100 900 560]);
subplot(2,1,1);
lag = (-12:12)';
plot(lag, max(dB(r/13), -80), 'ko-', 'LineWidth', 1.2); hold on;   % zeros drawn at -80 dB
lagm = (1:numel(ybest))' - (numel(ybest)+1)/2;
plot(lagm, max(dB(ybest), -80), 'r.-');
ylim([-80 3]); grid on; xlabel('lag (chips)'); ylabel('output (dB)');
legend('Barker-13 matched (-22.3 dB)', 'LS mismatched, L = 39', 'Location', 'northeast');
title('Barker-13: matched vs least-squares mismatched filter');
subplot(2,1,2);
tt = ((1:Nfft)' - Nfft/2 - 1)/fs*1e6;
plot(tt, dB(eA), 'Color', [0.6 0.6 0.6]); hold on;
plot(tt, dB(envs{1}), 'k', tt, dB(envs{2}), 'b', tt, dB(envs{3}), 'r');
xlim([-4 4]); ylim([-70 3]); grid on; xlabel('lag (\mus)'); ylabel('envelope (dB)');
legend('A alone', 'A+B static', 'A+B, d = 5 \mum', 'A+B, d = 20 \mum', 'Location', 'northeast');
title('Golay N = 16 on RF: motion between firings breaks the cancellation');
