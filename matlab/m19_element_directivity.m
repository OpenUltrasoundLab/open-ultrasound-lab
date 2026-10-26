%% m19_element_directivity.m  --  element factor, acceptance angle, F#min, and array pattern = element x array factor
% Purpose : (1) directivity of a strip element of width w: rigid baffle sinc(w sin(th)/lambda),
%           soft baffle sinc(...)*cos(th). Find the -6 dB half-angle by bisection and the
%           minimum useful F-number F#min = 1/(2 tan th6) for several widths and frequencies.
%           (2) steered array pattern P(th) = E(th) * AF(sin th - sin th0) for a linear array
%           (d = lambda) and a phased array (d = lambda/2): grating-lobe angles and levels.
% Change  : widths w_mm, frequencies f_MHz, N, pitch, steering angles.
% Expect  : w = 0.27 mm at 5 MHz: 43.5 deg rigid / 36.0 deg soft, F#min 0.53 / 0.69;
%           at 7.5 MHz 27.3 / 25.4 deg, F#min 0.97 / 1.05; at 3 MHz the rigid model never
%           reaches -6 dB, soft 48.9 deg (F#min 0.44).
%           Linear d = lambda, w = 0.9d, steered 20 deg: main-lobe loss -1.9 dB, grating lobe
%           at -41.1 deg, -6.3 dB re main. Phased d = lambda/2 at 45 deg: no grating lobe,
%           main-lobe loss -4.5 dB (w = 0.45 lambda, soft).
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
clear; close all; rng(1);

c = 1540;                                            % m/s
sincf = @(x) (abs(x) < 1e-12) + (abs(x) >= 1e-12) .* sin(pi*x) ./ (pi*x + (abs(x) < 1e-12));
Erig  = @(th, wl) abs(sincf(wl*sin(th)));            % rigid baffle, wl = w/lambda
Esoft = @(th, wl) abs(sincf(wl*sin(th))) .* cos(th); % soft baffle

%% (1) -6 dB acceptance angle and F#min
w_mm  = [0.27 0.27 0.27 0.154 0.30];
f_MHz = [3    5    7.5  5     2.5 ];
fprintf('  w (mm)  f (MHz)   w/lam   th6 rigid  th6 soft   F#min rigid  F#min soft\n');
for k = 1:numel(w_mm)
    lam = c/(f_MHz(k)*1e6)*1e3;  wl = w_mm(k)/lam;  th6 = [90 90];
    for b = 1:2
        if b == 1, E = @(t) Erig(t, wl); else, E = @(t) Esoft(t, wl); end
        if E(pi/2 - 1e-9) < 0.5
            lo = 1e-6; hi = pi/2 - 1e-9;
            for it = 1:60, mid = (lo + hi)/2; if E(mid) > 0.5, lo = mid; else, hi = mid; end, end
            th6(b) = (lo + hi)/2*180/pi;
        end
    end
    Fm = 1 ./ (2*tan(th6*pi/180));  Fm(th6 >= 89.99) = 0;
    fprintf('  %6.3f  %6.1f   %6.3f   %8.1f  %8.1f   %10.2f  %10.2f\n', w_mm(k), f_MHz(k), wl, th6(1), th6(2), Fm(1), Fm(2));
end

th = linspace(-pi/2, pi/2, 1801);
figure('Name', 'Element directivity');
cols = {'b', 'r', 'k'}; fl = [3 5 7.5];
for k = 1:3
    lam = c/(fl(k)*1e6)*1e3; wl = 0.27/lam;
    plot(th*180/pi, 20*log10(max(Erig(th, wl), 1e-3)), [cols{k} '--']); hold on;
    plot(th*180/pi, 20*log10(max(Esoft(th, wl), 1e-3)), cols{k}, 'LineWidth', 1.5);
end
plot([-90 90], [-6 -6], 'k:'); ylim([-40 0]); xlim([-90 90]); grid on;
xlabel('angle \theta (deg)'); ylabel('directivity (dB)');
title('w = 0.27 mm at 3 / 5 / 7.5 MHz (dashed rigid, solid soft)');

%% (2) Steered array pattern: element factor x array factor
N = 64;  u = sin(th);
AF = @(u, dl) abs(sum(exp(1i*2*pi*dl*((0:N-1).' - (N-1)/2)*u), 1)) / N;   % dl = d/lambda
cases = {'linear d = lambda, 20 deg', 1.0, 20; 'linear d = lambda, 0 deg', 1.0, 0; 'phased d = lambda/2, 45 deg', 0.5, 45};
figure('Name', 'Array pattern');
for k = 1:size(cases, 1)
    dl = cases{k, 2}; th0 = cases{k, 3}*pi/180; wl = 0.9*dl;
    P = AF(u - sin(th0), dl) .* Esoft(th, wl);
    main = Esoft(th0, wl);
    fprintf('\n%s: main-lobe loss E(th0) = %.2f dB\n', cases{k, 1}, 20*log10(main));
    found = 0;
    for m = [-3 -2 -1 1 2 3]
        s = sin(th0) + m/dl;
        if abs(s) <= 1
            tg = asin(s); found = 1;
            Eg = max(Esoft(tg, wl), 1e-5);
            fprintf('   grating lobe m = %+d at %.1f deg, %.1f dB re main (element factor %.1f dB; -100 = endfire null)\n', ...
                m, tg*180/pi, 20*log10(Eg/main), 20*log10(Eg));
        end
    end
    if ~found, fprintf('   no grating lobe in visible space\n'); end
    subplot(3, 1, k); plot(th*180/pi, 20*log10(max(P, 1e-4)), 'b'); hold on;
    plot(th*180/pi, 20*log10(max(Esoft(th, wl), 1e-4)), 'r--');
    ylim([-50 0]); xlim([-90 90]); grid on; ylabel('dB'); title(cases{k, 1});
end
xlabel('angle \theta (deg)');
