%% m34_gsdf.m  --  DICOM Grayscale Standard Display Function: L(j), j(L), JND count, calibration LUT
% Purpose : DICOM PS3.14 (chapter 7) defines the GSDF as a rational function of x = ln j,
%             log10 L(j) = (a + c x + e x^2 + g x^3 + m x^4) / (1 + b x + d x^2 + f x^3 + h x^4 + k x^5),
%           j = 1 ... 1023 (JND index), L from 0.05 to about 4000 cd/m^2, and an inverse polynomial
%             j(L) = A + B y + C y^2 + ... + I y^8,   y = log10 L.
%           A display realizes N_JND = j(L'max) - j(L'min) just-noticeable differences, where
%           L' = L + L_amb includes the diffusely reflected room light. A GSDF-calibrated display
%           spaces the P-values equally in j: N_JND/(2^b - 1) JNDs per step. A gamma-2.2 display
%           (same Lmin, Lmax) wastes its darkest steps (well under 1 JND), spaces the lower grays
%           up to 4 JNDs apart (visible banding) and squeezes the bright steps to about half the
%           GSDF spacing.
%           Last part: build the P-value -> DDL lookup table that calibrates a gamma-2.2 panel
%           to the GSDF, with an 8-bit and a 10-bit internal drive, and count lost levels.
% Change  : Lmin, Lmax, ambient list Lamb, panel gamma, bit depths.
% Expect  : L(1) = 0.050, L(1023) = 3993 cd/m^2; inverse agrees with the forward law to 0.1 JND.
%           0.5-400 cd/m^2: 626 JNDs dark room, 615 at 0.2 cd/m^2 ambient, 557 at 2 cd/m^2;
%           8-bit GSDF: 2.41 JND/step at 0.2 cd/m^2; gamma 2.2: 0.10 JND/step at the black end,
%           4.2 near P = 27 and 1.3 at the white end. Calibrating through an 8-bit drive leaves
%           213 distinct DDLs of 256; a 10-bit drive keeps all 256.
% Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.

%% GSDF coefficients (DICOM PS3.14, chapter 7)
a = -1.3011877;   b = -2.5840191e-2; c = 8.0242636e-2;  d = -1.0320229e-1; e = 1.3646699e-1;
f = 2.8745620e-2; g = -2.5468404e-2; h = -3.1978977e-3; k = 1.2992634e-4;  m = 1.3635334e-3;
A = 71.498068; B = 94.593053; C = 41.912053; D = 9.8247004; E = 0.28175407;
F = -1.1878455; G = -0.18014349; H = 0.14710899; I = -0.017046845;
Lj = @(j) 10.^((a + c*log(j) + e*log(j).^2 + g*log(j).^3 + m*log(j).^4) ./ ...
             (1 + b*log(j) + d*log(j).^2 + f*log(j).^3 + h*log(j).^4 + k*log(j).^5));
jL = @(L) A + B*log10(L) + C*log10(L).^2 + D*log10(L).^3 + E*log10(L).^4 + ...
          F*log10(L).^5 + G*log10(L).^6 + H*log10(L).^7 + I*log10(L).^8;

jj = 1:1023;
fprintf('L(1) = %.4f cd/m^2, L(1023) = %.1f cd/m^2\n', Lj(1), Lj(1023));
fprintf('inverse polynomial vs forward law: max |j(L(j)) - j| = %.3f JND\n', max(abs(jL(Lj(jj)) - jj)));
fprintf('luminance step per JND: %.2f %% at 1 cd/m^2, %.2f %% at 100 cd/m^2\n', ...
        100*(Lj(jL(1) + 1)/Lj(jL(1)) - 1), 100*(Lj(jL(100) + 1)/Lj(jL(100)) - 1));

%% JND count vs luminance range and ambient
Lmin = 0.5;  Lmax = 400;  Lamb = [0 0.2 2 10 50];
fprintf('\nLmin = %.2f, Lmax = %.0f cd/m^2\n  L_amb   LR''    N_JND   JND/step 8 b   10 b\n', Lmin, Lmax);
for q = 1:numel(Lamb)
    N = jL(Lmax + Lamb(q)) - jL(Lmin + Lamb(q));
    fprintf('  %5.1f  %6.1f  %6.1f   %6.3f      %6.3f\n', Lamb(q), (Lmax + Lamb(q))/(Lmin + Lamb(q)), N, N/255, N/1023);
end
N = jL(350) - jL(1);
fprintf('L''min = 1, L''max = 350 cd/m^2 (LR'' = 350): %.1f JNDs\n', N);

%% GSDF vs gamma 2.2: JNDs per gray step
gam = 2.2;  amb = [0.2 2];
figure('Position', [60 60 1100 380]);
subplot(1, 3, 1);
semilogy(jj, Lj(jj), 'k', 'LineWidth', 1.5); grid on;
xlabel('JND index j'); ylabel('L (cd/m^2)'); title('GSDF L(j)');
for q = 1:2
    for bits = [8 10]
        pmax = 2^bits - 1;  p = (0:pmax)/pmax;
        Lg = Lmin + (Lmax - Lmin)*p.^gam + amb(q);
        dJg = diff(jL(Lg));
        Ng = jL(Lmax + amb(q)) - jL(Lmin + amb(q));
        fprintf('ambient %.1f, %2d bit: GSDF %.3f JND/step; gamma 2.2 min %.3f max %.3f, steps > 1 JND: %d of %d\n', ...
                amb(q), bits, Ng/pmax, min(dJg), max(dJg), sum(dJg > 1), pmax);
        if bits == 8
            subplot(1, 3, 1 + q);
            plot(p(1:end-1), dJg, 'r', p([1 end-1]), [1 1]*Ng/pmax, 'b', p([1 end-1]), [1 1], 'k:', 'LineWidth', 1.4);
            grid on; xlabel('P-value / 255'); ylabel('JNDs per step');
            title(sprintf('8 bit, ambient %.1f cd/m^2', amb(q)));
            legend('gamma 2.2', 'GSDF', '1 JND', 'Location', 'northeast');
        end
    end
end

%% Calibration LUT: GSDF on a gamma-2.2 panel with an 8-bit or 10-bit drive
Lamb0 = 0.2;  jmin = jL(Lmin + Lamb0);  jmax = jL(Lmax + Lamb0);
P = 0:255;
Ltarget = Lj(jmin + P/255*(jmax - jmin)) - Lamb0;          % luminance the panel must emit
fprintf('\ncalibration LUT (P-value 0..255 -> DDL), panel gamma %.1f, ambient %.1f cd/m^2\n', gam, Lamb0);
for nb = [8 10]
    dmax = 2^nb - 1;
    DDL = round(dmax * ((Ltarget - Lmin)/(Lmax - Lmin)).^(1/gam));   % invert the panel characteristic
    Lout = Lmin + (Lmax - Lmin)*(DDL/dmax).^gam + Lamb0;
    dJ = diff(jL(Lout));
    err = 100*(dJ/((jmax - jmin)/255) - 1);
    fprintf('  %2d-bit drive: distinct DDLs %3d of 256, JND/step %.2f..%.2f (target %.2f), max step error %.0f %%\n', ...
            nb, numel(unique(DDL)), min(dJ), max(dJ), (jmax - jmin)/255, max(abs(err)));
end
