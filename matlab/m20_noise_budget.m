%% m20_noise_budget.m
% Element, cable and LNA: signal transfer, noise figure and noise match (Module 20.3-20.4).
%   Element: Butterworth-Van Dyke model (C0 in parallel with Rm, Lm, Cm),
%   resonant at f0 with effective coupling keff^2 and loaded Q.
%   Load: cable capacitance + LNA input capacitance, optional 200 ohm
%   termination (resistive or active) and optional shunt tuning inductor.
%   Noise at the LNA input: 4kT Re(Zs) + en^2 + in^2 |Zs|^2, with Zs the
%   impedance seen by the LNA. NF is referred to the element's own thermal
%   noise passed through the same divider as the signal.
% What to change: C0, keff2, Q, cable capacitance, en, in.
% Expected output: |Ze| = 637 ohm at -37 deg; cable 0/50/150/300 pF gives
%   |H| = -1.61/-6.28/-12.00/-16.83 dB and NF = 0.88/1.60/3.72/6.90 dB;
%   tuned: +0.63 dB, NF 1.28 dB; Ropt = 533 ohm, NFmin = 0.61 dB.
clear; close all; rng(1);

%% Parameters
kB  = 1.380649e-23;  T0 = 290;      % Boltzmann constant (J/K), noise-figure temperature (K)
f0  = 5e6;                          % element resonance (Hz)
C0  = 30e-12;                       % clamped capacitance (F)
k2  = 0.4;                          % keff^2 = Cm/(C0 + Cm)
Q   = 2;                            % loaded mechanical Q (broadband element)
Cin = 15e-12;                       % LNA input capacitance (F)
en  = 0.8e-9;                       % LNA voltage noise (V/sqrt(Hz))
in  = 1.5e-12;                      % LNA current noise (A/sqrt(Hz))
QL  = 30;                           % tuning-inductor Q
Rt  = 200;                          % termination resistance (ohm)

Cm = C0*k2/(1-k2); w0 = 2*pi*f0;
Rm = 1/(w0*Cm*Q);  Lm = 1/(w0^2*Cm);
fprintf('BVD: C0 = %.0f pF, Cm = %.1f pF, Lm = %.2f uH, Rm = %.1f ohm\n', C0*1e12, Cm*1e12, Lm*1e6, Rm);

%% Model as anonymous functions of frequency f (vector) and options
Ze   = @(f) 1 ./ (1./(Rm + 1i*2*pi*f*Lm + 1./(1i*2*pi*f*Cm)) + 1i*2*pi*f*C0);
% load admittance: cable + Cin, termination flag term (0 none, 1 resistive, 2 active), tuning flag tu
Yl   = @(f,Cc,term,tu) 1i*2*pi*f*(Cc+Cin) + (term>0)/Rt + ...
        tu ./ (w0*(1/(w0^2*(C0+Cc+Cin)))/QL + 1i*2*pi*f*(1/(w0^2*(C0+Cc+Cin))));
Hf   = @(f,Cc,term,tu) (1./Yl(f,Cc,term,tu)) ./ (Ze(f) + 1./Yl(f,Cc,term,tu));
Zsf  = @(f,Cc,term,tu) 1 ./ (1./Ze(f) + Yl(f,Cc,term,tu));
% total noise PSD at the LNA input (V^2/Hz); active termination: feedback resistor (1+A) = 11 times larger
Ntot = @(f,Cc,term,tu) 4*kB*T0*real(Zsf(f,Cc,term,tu)) + en^2 + in^2*abs(Zsf(f,Cc,term,tu)).^2 ...
        - (term==2)*(1-1/11)*4*kB*T0*abs(Zsf(f,Cc,term,tu)).^2/Rt;
Nel  = @(f,Cc,term,tu) 4*kB*T0*real(Ze(f)).*abs(Hf(f,Cc,term,tu)).^2;
NF   = @(f,Cc,term,tu) 10*log10(Ntot(f,Cc,term,tu)./Nel(f,Cc,term,tu));

%% Spot values at f0
z = Ze(f0);
fprintf('Element at f0: |Ze| = %.0f ohm, angle %.0f deg, Re = %.0f ohm, thermal noise %.2f nV/rtHz\n', ...
  abs(z), angle(z)*180/pi, real(z), sqrt(4*kB*T0*real(z))*1e9);
fprintf('\n%-34s %9s %9s %12s\n', 'configuration', '|H| (dB)', 'NF (dB)', 'noise nV/rtHz');
cases = {0,0,0,'no cable, high-Z'; 50e-12,0,0,'50 pF cable, high-Z'; 150e-12,0,0,'150 pF cable, high-Z'; ...
         300e-12,0,0,'300 pF cable, high-Z'; 150e-12,0,1,'150 pF, tuned (QL = 30)'; ...
         150e-12,1,0,'150 pF, 200 ohm resistor'; 150e-12,2,0,'150 pF, 200 ohm active (x11)'};
for k = 1:size(cases,1)
  Cc = cases{k,1}; term = cases{k,2}; tu = cases{k,3};
  fprintf('%-34s %9.2f %9.2f %12.2f\n', cases{k,4}, 20*log10(abs(Hf(f0,Cc,term,tu))), NF(f0,Cc,term,tu), sqrt(Ntot(f0,Cc,term,tu))*1e9);
end
fprintf('capacitive-divider estimate for 150 pF: %.2f dB\n', 20*log10(C0/(C0+150e-12+Cin)));
fprintf('tuning inductor for 150 pF: L = %.2f uH\n', 1/(w0^2*(C0+150e-12+Cin))*1e6);

%% Noise match for a resistive source
Ropt  = en/in;  NFmin = 10*log10(1 + en*in/(2*kB*T0));
fprintf('\nResistive source: Ropt = en/in = %.0f ohm, NFmin = %.2f dB\n', Ropt, NFmin);
fprintf('sqrt(4kTR): 50 ohm %.3f, 1 kohm %.3f nV/rtHz; en equals a %.1f ohm resistor\n', ...
  sqrt(4*kB*T0*50)*1e9, sqrt(4*kB*T0*1000)*1e9, en^2/(4*kB*T0));
% cascade with a second stage of 5 nV/rtHz
for A1dB = [6 20]
  fprintf('LNA gain %2d dB: input-referred total %.2f nV/rtHz (e2 = 5 nV/rtHz)\n', A1dB, sqrt(en^2 + (5e-9)^2/10^(A1dB/10))*1e9);
end

%% Figures
f = linspace(0.5e6, 10e6, 400);
figure('Name','m23 receive transfer and NF');
subplot(2,1,1);
plot(f/1e6, 20*log10(abs(Hf(f,0,0,0))), f/1e6, 20*log10(abs(Hf(f,150e-12,0,0))), f/1e6, 20*log10(abs(Hf(f,300e-12,0,0))), f/1e6, 20*log10(abs(Hf(f,150e-12,0,1))), '--');
grid on; ylabel('|H| (dB)'); title('V_{LNA}/V_{oc}: element + cable + LNA input'); legend('0 pF','150 pF','300 pF','150 pF tuned','Location','southwest'); ylim([-40 10]);
subplot(2,1,2);
plot(f/1e6, NF(f,0,0,0), f/1e6, NF(f,150e-12,0,0), f/1e6, NF(f,300e-12,0,0), f/1e6, NF(f,150e-12,0,1), '--');
grid on; ylim([0 20]); xlabel('frequency (MHz)'); ylabel('NF (dB)'); title('noise figure');
figure('Name','m23 noise match');
R = logspace(1, 4, 300);
semilogx(R, 10*log10(1 + (en^2 + in^2*R.^2)./(4*kB*T0*R)), Ropt, NFmin, 'o'); grid on;
xlabel('source resistance R_s (\Omega)'); ylabel('NF (dB)'); title(sprintf('NF vs R_s, R_{opt} = %.0f \\Omega', Ropt)); ylim([0 10]);
