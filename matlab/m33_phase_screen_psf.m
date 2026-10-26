%% m33_phase_screen_psf.m
% Thin near-field phase screen: Strehl ratio, PSF broadening and the loss of
% spatial coherence, by Monte Carlo against the closed forms.
% Purpose : A Gaussian random delay screen tau(x) (rms sigma_tau, correlation
%           rho(dx) = exp(-dx^2/l^2)) sits on a 128-element array (pitch
%           0.3 mm, f0 = 5 MHz). Focal-plane (Fraunhofer) PSF of the focused
%           aperture, one-way p(x) = sum_i exp(j phi_i) exp(-j k x_i x / z).
%           (A) On-axis Strehl S = |sum exp(j phi_i)|^2 / N^2 against
%                exp(-sigma_phi^2) (infinite aperture) and the exact mean
%                E[S] = (1/N^2) sum_ij exp(-sigma_phi^2 (1 - rho_ij)).
%           (B) Mean PSF intensity: Monte Carlo against the exact
%                E|p(x)|^2 = sum_ij exp(-j k x (x_i - x_j)/z) exp(-sigma^2(1-rho_ij)),
%                and its -6 dB width (intensity 1/4) against the halo of the
%                large-sigma limit, a Gaussian of rms s = sqrt(2) c sigma_tau z / l
%                (-6 dB width 3.33 s).
%           (C) Receive coherence of speckle (ideal transmit, screen on
%                receive): R(m) = (1 - m/N) exp(-sigma_phi^2 (1 - rho(m p))).
% Change  : sig_list (ns), ell (mm), f0, z.
% Expect  : Strehl follows exp(-sigma_phi^2) for small sigma and then the
%           slower sqrt(pi) l/(sigma_phi D); the mean PSF keeps its core until
%           sigma_phi ~ 1 rad and then becomes the halo; coherence falls as predicted.
% Runtime : about 5 s in Octave.
clear; close all; rng(1);

c = 1540; f0 = 5e6; lam = c/f0; k = 2*pi/lam;
N = 128; p = 0.3e-3; xe = ((0:N-1) - (N-1)/2)*p; D = N*p;   % element positions (m), aperture 38.4 mm
z = 40e-3; ell = 4e-3;                                       % focal depth, screen correlation length (m)
sig_list = [10 20 30 40 60 80]*1e-9;                         % rms delay (s)
R = 400;                                                     % realizations

% screen generator: white noise on a fine grid convolved with a Gaussian kernel
dxs = p/4; xs = (-60e-3:dxs:60e-3);                          % fine grid (m), much longer than the aperture
g = exp(-2*(((-3*ell):dxs:(3*ell))/ell).^2); g = g/sqrt(sum(g.^2));   % autocorrelation exp(-dx^2/l^2), unit variance
idx = round((xe - xs(1))/dxs) + 1;                           % element samples on the fine grid
screen = @(sig) sig*conv(randn(1, numel(xs)), g, 'same');    % one realization tau(xs), rms sig

DX = xe' - xe; RHO = exp(-(DX/ell).^2);                      % rho_ij
fprintf('Aperture %.1f mm, z = %.0f mm (F/%.2f), l = %.0f mm, lambda = %.3f mm\n', D*1e3, z*1e3, z/D, ell*1e3, lam*1e3);
fprintf('Large-sigma asymptote of E[S]: sqrt(pi) l / (sigma_phi D) (one correlation cell of N)\n\n');
fprintf(' sig_tau  sig_phi   S MC     E[S] exact  exp(-s^2)  asympt.   S2way MC  E[S]^2\n');
xf = linspace(-4e-3, 4e-3, 801);                             % lateral focal coordinate (m)
E = exp(-1j*k*xf'*xe/z);                                     % Fraunhofer kernel (801 x N)
w6 = zeros(size(sig_list)); w6x = w6; blur = w6; S_mc = w6; S_ex = w6;
for a = 1:numel(sig_list)
    sig = sig_list(a); sp = 2*pi*f0*sig;
    S = zeros(R, 1); S2 = S; I = zeros(numel(xf), 1);
    chk = 0;
    for r = 1:R
        t = screen(sig); ph = 2*pi*f0*t(idx);
        chk = chk + mean(t(idx).^2);
        S(r) = abs(sum(exp(1j*ph)))^2/N^2; S2(r) = S(r)^2;
        I = I + abs(E*exp(1j*ph(:))).^2/N^2;
    end
    I = I/R; Sex = sum(sum(exp(-sp^2*(1 - RHO))))/N^2;
    S_mc(a) = mean(S); S_ex(a) = Sex;
    fprintf('%5.0f ns  %5.3f    %.3f    %.3f       %.3f      %.3f     %.3f     %.3f\n', sig*1e9, sp, mean(S), Sex, exp(-sp^2), sqrt(pi)*ell/(sp*D), mean(S2), Sex^2);
    % exact mean PSF intensity
    Iex = real(sum((E*exp(-sp^2*(1 - RHO))).*conj(E), 2))/N^2;
    above = xf(I >= max(I)/4); w6(a) = max(above) - min(above);           % -6 dB (intensity 1/4) width, MC
    above = xf(Iex >= max(Iex)/4); w6x(a) = max(above) - min(above);      % exact
    blur(a) = sqrt(2)*c*sig*z/ell;                               % rms width s of the Gaussian intensity halo
    if a == 1, fprintf('   (rms of generated screen at the elements: %.2f ns)\n', sqrt(chk/R)*1e9); end
    if a == 4, I40 = I; Iex40 = Iex; end
end
Pid = abs(E*ones(N, 1)).^2/N^2; above = xf(Pid >= 1/4); w0 = max(above) - min(above);
fprintf('\n-6 dB width of the mean one-way PSF intensity (mm); unaberrated %.3f mm (1.21 lambda F = %.3f)\n', w0*1e3, 1.207*lam*z/D*1e3);
fprintf(' sig_tau   MC     exact   halo 3.33 s, s = sqrt(2) c sig z / l (-6 dB width of the Gaussian halo)\n');
for a = 1:numel(sig_list)
    fprintf('%5.0f ns  %.3f   %.3f   %.3f\n', sig_list(a)*1e9, w6(a)*1e3, w6x(a)*1e3, 3.330*blur(a)*1e3);
end

% (C) receive coherence of speckle with an ideal transmit and the screen on receive
sig = 40e-9; sp = 2*pi*f0*sig; Ns = 1500; Rc = 300; lags = 0:40; Rm = zeros(size(lags)); Pa = Rm; Pb = Rm;
for r = 1:Rc
    t = screen(sig); ph = 2*pi*f0*t(idx);
    xsc = (rand(Ns, 1) - 0.5)*30*lam*z/D; b = sin(pi*xsc*D/(lam*z) + eps)./(pi*xsc*D/(lam*z) + eps);   % ideal transmit beam amplitude
    a_s = (randn(Ns, 1) + 1j*randn(Ns, 1)).*b;
    s = exp(1j*ph(:)).*(exp(-1j*k*xe'*xsc'/z)*a_s);          % received aperture vector
    for m = lags
        Rm(m+1) = Rm(m+1) + real(sum(s(1:N-m).*conj(s(1+m:N))));   % accumulate numerator over realizations
        Pa(m+1) = Pa(m+1) + sum(abs(s(1:N-m)).^2); Pb(m+1) = Pb(m+1) + sum(abs(s(1+m:N)).^2);
    end
end
Rm = Rm./sqrt(Pa.*Pb);                                        % ensemble (depth-averaged) coherence estimate
Rth = (1 - lags/N).*exp(-sp^2*(1 - exp(-(lags*p/ell).^2)));
R0 = 1 - lags/N;
fprintf('\nReceive coherence, 40 ns screen (sigma_phi = %.2f rad), l = 4 mm\n lag   MC     formula   no screen\n', sp);
for m = [1 2 5 10 13 20 40]
    fprintf('%3d   %.3f   %.3f     %.3f\n', m, Rm(m+1), Rth(m+1), R0(m+1));
end

figure('Position', [50 50 1000 330]);
subplot(1, 3, 1); sl = linspace(0, 90, 200)*1e-9; spl = 2*pi*f0*sl; Sl = zeros(size(sl));
for q = 1:numel(sl), Sl(q) = sum(sum(exp(-spl(q)^2*(1 - RHO))))/N^2; end
plot(sl*1e9, exp(-spl.^2), 'k--', sl*1e9, Sl, 'b-', sig_list*1e9, S_mc, 'ro'); xlabel('\sigma_\tau (ns)'); ylabel('Strehl S'); title('On-axis Strehl, one way'); legend('exp(-\sigma_\phi^2)', 'exact mean', 'Monte Carlo');
subplot(1, 3, 2); plot(xf*1e3, 10*log10(Pid), 'k', xf*1e3, 10*log10(I40), 'r', xf*1e3, 10*log10(Iex40), 'b--'); ylim([-40 0]); xlabel('x (mm)'); ylabel('dB'); title('Mean PSF intensity, 40 ns'); legend('no screen', 'MC', 'exact');
subplot(1, 3, 3); plot(lags, R0, 'k', lags, Rth, 'b-', lags, Rm, 'ro'); xlabel('lag m (elements)'); ylabel('coherence'); title('Receive coherence, 40 ns'); legend('no screen', 'formula', 'MC');
