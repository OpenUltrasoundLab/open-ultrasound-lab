%% m19_rca_beam.m  --  row-column addressed array: separable PSF and edge-wave ghosts
% Purpose : (1) focal-plane (Fourier) PSF of an N x N row-column array (transmit on rows: focus
%           in y only; receive on columns: focus in x only) versus a fully addressed matrix
%           focused on transmit and receive: PSF_RCA = R(x) T(y), PSF_2D = Tx Rx (x) . Ty Ry (y).
%           (2) pulse-echo signal of one transmit row and one receive column for one scatterer:
%           s = p * h_row * h_col, with h the spatial impulse response of a strip (1/r line of
%           point sources). Uniform strips have steps at their ends (edge waves -> ghost echoes);
%           roll-off strips (Hann taper added beyond each end, as in Rasmussen et al. 2015)
%           remove them.
% Change  : apodization, N, scatterer position (xs, ys, zs), roll-off width.
% Expect  : uniform: axis PSL RCA -13.3 dB vs full -26.5 dB; diagonal -26.5 vs -53.0 dB;
%           -6 dB width 1.21 vs 0.89 lambda F#. Hann: axis -31.5 vs -62.9 dB.
%           Echo at (8, 3, 30) mm: uniform-strip edge echoes about -22, -27, -30 dB after the main
%           echo; with a 16 lambda roll-off the highest is below about -43 dB. Centred scatterer
%           (0, 0, 30): edge echo -13.3 dB (both ends arrive together).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

%% (1) Focal-plane PSF
N = 64; x = linspace(-4, 4, 1601);                  % lateral position in units of lambda F#
n = (0:N-1) - (N-1)/2;
apods = {'uniform', ones(1, N); 'Hann', 0.5 - 0.5*cos(2*pi*((0:N-1) + 0.5)/N)};
fprintf('(1) N = %d per side: RCA uses %d channels, fully addressed %d\n', N, 2*N, N^2);
for a = 1:2
    w = apods{a, 2};
    p = abs(exp(1i*2*pi*x(:)*n/N) * w(:)).' / sum(w);    % one-way focal pattern
    xx = x(x >= 0); pp = p(x >= 0); k = 2; while pp(k+1) < pp(k), k = k + 1; end
    psl = max(pp(k:end));
    w1 = 2*xx(find(pp < 0.5, 1)); w2 = 2*xx(find(pp.^2 < 0.5, 1));
    fprintf('   %-8s axis PSL RCA %6.1f dB, full %6.1f dB | diagonal RCA %6.1f, full %6.1f dB | -6 dB width RCA %.2f, full %.2f lambda F#\n', ...
        apods{a, 1}, 20*log10(psl), 40*log10(psl), 40*log10(psl), 80*log10(psl), w1, w2);
    if a == 1, p1 = p; end
end
xs2 = x(1:8:end); q = p1(1:8:end);
figure('Name', 'RCA PSF');
subplot(2, 2, 1); imagesc(xs2, xs2, max(20*log10(q.' * q), -50)); axis image; caxis([-50 0]); title('RCA: R(x) T(y)'); xlabel('x (\lambda F#)'); ylabel('y');
subplot(2, 2, 2); imagesc(xs2, xs2, max(20*log10((q.^2).' * q.^2), -50)); axis image; caxis([-50 0]); title('fully addressed'); xlabel('x (\lambda F#)'); colorbar;
subplot(2, 1, 2); plot(x, 20*log10(p1 + 1e-9), 'r', x, 40*log10(p1 + 1e-9), 'b'); ylim([-60 0]); grid on;
xlabel('x (\lambda F#)'); ylabel('dB'); legend('RCA axis (one-way)', 'full axis = RCA diagonal');

%% (2) Edge waves: one row (along x at y = 0) transmits, one column (along y at x = 0) receives
c = 1540e3; f0 = 5e6; lam = c/f0;                    % mm/s, Hz, mm
L0 = 128*lam/2; fs = 250e6; dt = 1/fs; nx = 4000;   % strip length (mm), sampling
sig = 1/(2*pi*0.6*f0/2.3548);                        % 60 % (-6 dB) bandwidth Gaussian pulse
tp = (-200:200)*dt; pulse = exp(-tp.^2/(2*sig^2)) .* cos(2*pi*f0*tp);
Pn = 8192; pre = 250;
cases = [8 3 30; 0 0 30];
figure('Name', 'RCA edge waves');
for ci = 1:2
    sc = cases(ci, :);
    for rw = [0 16*lam]
        L = L0 + 2*rw; s = -L/2 + ((0:nx-1) + 0.5)*L/nx; e = L/2 - abs(s);
        a = ones(1, nx); if rw > 0, t = e < rw; a(t) = 0.5 - 0.5*cos(pi*e(t)/rw); end
        rt = sqrt((s - sc(1)).^2 + sc(2)^2 + sc(3)^2);  rr = sqrt((s - sc(2)).^2 + sc(1)^2 + sc(3)^2);
        wt = a*(L/nx)./rt; wr = a*(L/nx)./rr;
        it = round((rt - min(rt))/c/dt) + pre + 1; ir = round((rr - min(rr))/c/dt) + 1;
        ht = accumarray(it(:), wt(:), [Pn 1]).'; hr = accumarray(ir(:), wr(:), [Pn 1]).';
        pk = zeros(1, Pn); pk(1:201) = pulse(201:401); pk(end-199:end) = pulse(1:200);   % centred, circular
        S = fft(ht) .* fft(hr) .* fft(pk);
        H = zeros(1, Pn); H(1) = 1; H(2:Pn/2) = 2; H(Pn/2+1) = 1;
        env = abs(ifft(S .* H)); env = env / max(env);
        tt = ((1:Pn) - pre - 1)*dt*1e6;              % us after the closest-point arrival
        % predicted edge arrivals (us): combinations of closest point and the two physical ends
        rtE = [min(rt), sqrt((-L0/2 - sc(1))^2 + sc(2)^2 + sc(3)^2), sqrt((L0/2 - sc(1))^2 + sc(2)^2 + sc(3)^2)];
        rrE = [min(rr), sqrt((-L0/2 - sc(2))^2 + sc(1)^2 + sc(3)^2), sqrt((L0/2 - sc(2))^2 + sc(1)^2 + sc(3)^2)];
        [A1, B1] = meshgrid(rtE - rtE(1), rrE - rrE(1)); tc = (A1(:) + B1(:)).'/c*1e6; tc = tc(tc > 0.45);
        best = -Inf;
        for q = tc, kk = abs(tt - q) <= 0.1; best = max(best, 20*log10(max(env(kk)))); end
        if rw == 0, lab = 'uniform'; else, lab = sprintf('roll-off %g lambda', rw/lam); end
        fprintf('(2) scatterer (%g, %g, %g) mm, %-16s: highest edge echo after 0.45 us %6.1f dB\n', sc(1), sc(2), sc(3), lab, best);
        subplot(2, 1, ci); plot(tt, 20*log10(env + 1e-9)); hold on;
    end
    xlim([-1 4]); ylim([-70 0]); grid on; ylabel('envelope (dB)');
    title(sprintf('one row Tx, one column Rx, scatterer (%g, %g, %g) mm', sc)); legend('uniform', 'roll-off 16\lambda');
end
xlabel('time after main echo (\mus)');
