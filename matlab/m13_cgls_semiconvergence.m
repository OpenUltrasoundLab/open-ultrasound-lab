%% m13_cgls_semiconvergence.m
% Iteration count as the regularization parameter: CGLS and Landweber.
% Purpose : Same RF deconvolution problem as m13_tikhonov_deconv_lcurve.m.
%           Run CGLS (conjugate gradients on the normal equations, never
%           forming A'A) and Landweber x <- x + w A'(b - A x), w = 1/sigma_max^2.
%           The error against the band-limited truth first falls, then rises
%           (semiconvergence). Landweber after k steps has filter factors
%           1 - (1 - w sigma^2)^k, half-power at sigma^2 ~ ln2/k, so it acts
%           like Tikhonov with lambda^2 ~ ln2/k.
% Change  : snr_dB; kmax; the discrepancy factor tau.
% Expect  : at 25 dB the best CGLS iterate is reached in tens of iterations,
%           Landweber needs about ln2/lambda_opt^2 (over a thousand), and
%           both reach about the oracle-Tikhonov error.
% Runtime : about 1 s.
clear; close all; rng(1);

N = 512; fs = 40e6; f0 = 5e6; bw = 0.6;
sf = bw*f0/2/sqrt(2*log(2)); st = 1/(2*pi*sf);
t = ((0:N-1)' - N/2)/fs;
h = circshift(exp(-t.^2/(2*st^2)).*cos(2*pi*f0*t), -N/2);
H = fft(h); H = H/max(abs(H));
Aop  = @(v) real(ifft(H.*fft(v)));
Atop = @(v) real(ifft(conj(H).*fft(v)));        % adjoint (correlation with the pulse)

x = 0.08*randn(N,1);
pos = [70 150 160 260 330 345 420]; amp = [1 -0.7 0.6 0.9 -0.5 0.5 0.8];
x(pos) = x(pos) + amp(:);
band = abs(H) >= 0.01; X = fft(x); X(~band) = 0; xB = real(ifft(X));
b0 = Aop(x);
w = randn(N,1);   % one draw reused so the three SNRs differ only in level

tau = 1.0; kmax = 400;
for snr_dB = [40 25 10]
    n = w*norm(b0)/norm(w)*10^(-snr_dB/20); sig = norm(n)/sqrt(N); b = b0 + n;
    % ---- CGLS
    xk = zeros(N,1); r = b; s = Atop(r); p = s; g = s'*s;
    Ec = zeros(kmax,1); Rc = Ec;
    for k = 1:kmax
        q = Aop(p); a = g/(q'*q);
        xk = xk + a*p; r = r - a*q;
        s = Atop(r); gn = s'*s; p = s + (gn/g)*p; g = gn;
        Ec(k) = norm(xk - xB)/norm(xB); Rc(k) = norm(r);
    end
    [ec, kc] = min(Ec); kdp = find(Rc <= tau*sqrt(N)*sig, 1);
    % ---- Landweber, closed form in the Fourier domain (w = 1)
    ks = unique(round(logspace(0, 5, 300))); El = zeros(size(ks)); B = fft(b);
    Hs = H; Hs(abs(Hs) < 1e-300) = 1;
    for q = 1:numel(ks)
        F = 1 - (1 - abs(H).^2).^ks(q);
        El(q) = norm(real(ifft(F./Hs.*B)) - xB)/norm(xB);
    end
    [el, ql] = min(El);
    % ---- Tikhonov oracle for comparison
    lams = logspace(-4, 0.5, 300); Et = zeros(size(lams));
    for q = 1:numel(lams)
        Et(q) = norm(real(ifft(conj(H)./(abs(H).^2 + lams(q)^2).*B)) - xB)/norm(xB);
    end
    [et, qt] = min(Et);
    fprintf('SNR %2d dB: CGLS best k = %3d (err %.3f), discrepancy stop k = %3d (err %.3f), err at k=%d: %.2f\n', ...
        snr_dB, kc, ec, kdp, Ec(kdp), kmax, Ec(end));
    fprintf('           Landweber best k = %5d (err %.3f); Tikhonov oracle lambda = %.4f (err %.3f); ln2/lambda^2 = %.0f\n', ...
        ks(ql), el, lams(qt), et, log(2)/lams(qt)^2);
    if snr_dB == 25, Ec25 = Ec; El25 = El; et25 = et; kc25 = kc; kdp25 = kdp; end
end

figure('Position', [50 50 900 380]);
subplot(1,2,1);
semilogx(1:kmax, Ec25, 'b-'); hold on; plot(kc25, Ec25(kc25), 'bo'); plot(kdp25, Ec25(kdp25), 'ms');
plot([1 kmax], [et25 et25], 'k:'); ylim([0 1.5]);
xlabel('CGLS iteration k'); ylabel('relative error'); title('CGLS, 25 dB: semiconvergence');
legend('error', 'best', 'discrepancy stop', 'Tikhonov oracle', 'Location', 'northwest');
subplot(1,2,2);
semilogx(ks, El25, 'r-'); hold on; plot([1 1e5], [et25 et25], 'k:'); ylim([0 1.5]);
xlabel('Landweber iteration k'); ylabel('relative error'); title('Landweber (\omega = 1/\sigma_{max}^2), 25 dB');
