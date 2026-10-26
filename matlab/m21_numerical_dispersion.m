%% m21_numerical_dispersion.m
% 1D time stepping with three schemes over 100 wavelengths, compared with their dispersion relations.
%   staggered leapfrog:  sin(w dt/2) = CFL * S(k dx)
%   FD2  S = sin(q/2);  FD4  S = 9/8 sin(q/2) - 1/24 sin(3q/2);  k-space PSTD: exact (w = c k)
% What it shows: the measured phase lag at f0 after 100 lambda equals the dispersion-relation
%   prediction; FD2 at 8 PPW is more than two periods late, FD4 about 0.06 cycle (~20 deg) early, k-space exact.
% What to change: PPW, CFL, distance.
% Expected output (8 PPW, CFL 0.3): FD2 -2.34 cycles, FD4 +0.06, k-space 0.00 (predicted = measured).
% Units: wavelength lambda0 = 1, c = 1, rho = 1.  Base MATLAB / Octave.
% These are small 1D solvers written for this course, not k-Wave; the k-space scheme uses the
%   same correction idea (Tabei, Mast & Waag 2002) on a periodic 1D grid.

clear; close all; rng(1);
PPW = 8; CFL = 0.3; Ldist = 100;            % points per wavelength at f0, Courant number, distance (lambda)
Ndom = 128;                                 % domain length in wavelengths (periodic)
N = Ndom*PPW; dx = 1/PPW; dt = CFL*dx; nsteps = round(Ldist/dt); T = nsteps*dt;
x  = (0:N-1)*dx;                            % pressure nodes
xv = x + dx/2;                              % velocity nodes (staggered)
f0 = 1; B = 0.6; sf = B*f0/(2*sqrt(2*log(2))); st = 1/(2*pi*sf);
x0 = 20;                                    % initial pulse centre
pulse = @(xx) exp(-(xx - x0).^2/(2*st^2)) .* sin(2*pi*f0*(xx - x0));
k  = 2*pi*[0:N/2, -N/2+1:-1]/(N*dx);        % FFT wavenumbers
sinc_ = @(z) (abs(z) < 1e-12) + (abs(z) >= 1e-12) .* sin(z) ./ (z + (abs(z) < 1e-12));
q0 = 2*pi/PPW; k0bin = Ndom + 1;            % f0 sits exactly on bin Ndom (domain = 128 lambda)

schemes = {'FD2', 'FD4', 'k-space'};
Sfun = {@(q) sin(q/2), @(q) 9/8*sin(q/2) - 1/24*sin(3*q/2), []};
figure(1); clf;
for is = 1:3
    p = pulse(x);
    v = pulse(xv + dt/2);                   % right-going p(x - ct): v = p/(rho c) at t = -dt/2
    tic;
    switch is
        case 1
            for n = 1:nsteps
                v = v - CFL*(circshift(p, [0 -1]) - p);
                p = p - CFL*(v - circshift(v, [0 1]));
            end
        case 2
            A1 = 9/8; A2 = -1/24;
            for n = 1:nsteps
                v = v - CFL*(A1*(circshift(p, [0 -1]) - p) + A2*(circshift(p, [0 -2]) - circshift(p, [0 1])));
                p = p - CFL*(A1*(v - circshift(v, [0 1])) + A2*(circshift(v, [0 -1]) - circshift(v, [0 2])));
            end
        case 3
            kap = sinc_(k*dt/2);            % k-space correction, c_ref = 1
            sh = exp(1i*k*dx/2);            % staggering shift
            for n = 1:nsteps
                v = v - dt*real(ifft(1i*k.*kap.*sh .* fft(p)));
                p = p - dt*real(ifft(1i*k.*kap.*conj(sh) .* fft(v)));
            end
    end
    el = toc;
    % predicted phase error at f0 (cycles): (w_num - c k0) T / 2pi
    if is < 3
        wn = 2*asin(CFL*Sfun{is}(q0))/dt;
    else
        wn = 2*pi*f0;
    end
    pred = (wn - 2*pi*f0)*T/(2*pi);
    % measured: phase of bin f0 relative to the exact translation
    P0 = fft(pulse(x)); P1 = fft(p);
    expct = P0(k0bin) * exp(-1i*2*pi*f0*T);
    meas = angle(P1(k0bin)/expct)/(2*pi);   % cycles, wrapped to (-0.5, 0.5]
    wrapP = pred - round(pred);
    fprintf('%-8s predicted %+7.3f cycles (wrapped %+6.3f), measured (wrapped) %+6.3f, travel-time error %+6.3f %%, %.1f s\n', ...
        schemes{is}, pred, wrapP, -meas, -100*pred/Ldist, el);
    xs = mod(x - x0 - T + Ndom/2, Ndom) - Ndom/2;     % moving frame, centred
    [xs, io] = sort(xs);
    subplot(3,1,is); plot(xs, p(io), 'LineWidth', 1.2); hold on; plot(xs, pulse(xs + x0), 'k--');
    xlim([-8 4]); grid on; ylabel('p'); title(sprintf('%s after %d \\lambda at %g PPW, CFL %.2f', schemes{is}, Ldist, PPW, CFL));
end
xlabel('(x - ct)/\lambda_0');

% dispersion curves
ppw = linspace(3, 20, 200); qq = 2*pi./ppw;
figure(2); clf;
plot(ppw, 100*(2*asin(CFL*sin(qq/2))./(CFL*qq) - 1), 'r', ppw, 100*(2*asin(CFL*(9/8*sin(qq/2) - 1/24*sin(3*qq/2)))./(CFL*qq) - 1), 'Color', [0.7 0.45 0], 'LineWidth', 1.6); hold on;
plot(ppw, 100*(2*asin(min(1, CFL*qq/2))./(CFL*qq) - 1), 'b', ppw, 0*ppw, 'g--', 'LineWidth', 1.6);
ylim([-8 3]); grid on; xlabel('points per wavelength'); ylabel('phase-speed error (%)');
legend('FD2', 'FD4', 'pseudospectral leapfrog', 'k-space corrected', 'Location', 'southeast'); title(sprintf('Numerical phase speed, CFL %.2f', CFL));
fprintf('\nFD2 at 10 PPW: %.3f %%  (approx -(k dx)^2 (1 - CFL^2)/24 = %.3f %%)\n', 100*(2*asin(CFL*sin(pi/10))/(CFL*2*pi/10) - 1), -100*(2*pi/10)^2*(1 - CFL^2)/24);
