%% m20_pulser_spectra.m
% Pulser waveforms and their harmonic content (Module 20.2).
%   Builds unipolar, bipolar, 3-level (120 deg conduction), 5-level
%   (switching at 12 and 48 deg) and linear (Hann-windowed sine) bursts,
%   then measures the level at n*f0 by an exact DFT of the burst.
%   Part 1: ideal edges, 20-cycle bursts (steady-state harmonic content).
%   Part 2: slew-limited edges, I/C = 2 A / 300 pF at +-50 V.
%   Part 3: duty error D = 1/2 + Delta/T, simulated 2f0 versus sin(pi*Delta*f0).
% What to change: f0, V, Cload, Ipk, the conduction angle beta, the cycles.
% Expected output: bipolar 3f0 -9.5 dBc, 5f0 -14.0 dBc; 3-level no 3f0;
%   5-level no 3f0 or 5f0; duty error of 1.1 ns at 5 MHz: about -35 dBc,
%   matching theory to within 0.1 dB.
clear; close all; rng(1);

%% Parameters
f0    = 5e6;              % centre frequency (Hz)
V     = 50;               % rail voltage (V), drive is +-V
dt    = 1/(720*f0);       % 0.5 deg of phase per sample (0.278 ns)
N     = 2^15;             % samples (9.1 us window)
t     = (0:N-1)*dt;       % time axis (s)
t0    = 0.2e-6;           % burst start (s)
beta  = 120;              % 3-level conduction angle per half cycle (deg)
Cload = 300e-12;          % element + cable capacitance (F)
Ipk   = 2;                % pulser peak current (A)
slew  = Ipk/Cload;        % slew rate (V/s)

%% Level functions of the normalized phase u in [0,1)
lvUni  = @(u) double(u < 0.5);
lvBi   = @(u) 1 - 2*double(u >= 0.5);
lvTri  = @(u) double(abs(u-0.25) < beta/720) - double(abs(u-0.75) < beta/720);
lv5h   = @(a) 0.5*((a>=12)&(a<48)) + ((a>=48)&(a<132)) + 0.5*((a>=132)&(a<168));
lvFive = @(u) (u<0.5).*lv5h(u*360) - (u>=0.5).*lv5h(u*360-180);
% duty error e: the positive half cycle lasts (1/2 + e) of the period
dutyMap = @(fr,e) (fr < 0.5+e).*(fr/(0.5+e)*0.5) + (fr >= 0.5+e).*(0.5 + (fr-0.5-e)/(0.5-e)*0.5);
dft1    = @(x,f) abs(sum(x .* exp(-1i*2*pi*f*t)));   % exact DFT at one frequency
names   = {'unipolar', 'bipolar', '3-level', '5-level', 'linear'};

%% Waveform builder (ideal edges); type index k, cycles nc, duty error e
X = cell(5,1);
ph  = (t - t0)*f0;
for k = 1:5
  for pass = 1:2
    if pass == 1, nc = 20; else, nc = 4; end
    inb = (ph >= 0) & (ph < nc);
    fr  = ph - floor(ph);
    u   = dutyMap(fr, 0);
    switch k
      case 1, lv = lvUni(u);
      case 2, lv = lvBi(u);
      case 3, lv = lvTri(u);
      case 4, lv = lvFive(u);
      case 5, lv = sin(2*pi*ph) .* 0.5.*(1 - cos(2*pi*ph/nc));
    end
    x = V * lv .* inb;
    if pass == 1, X{k} = x; else, Xs{k} = x; end %#ok<SAGROW>
  end
end

%% Part 1: ideal edges, 20-cycle bursts
fprintf('Part 1: ideal edges, 20 cycles at %.1f MHz, +-%d V\n', f0/1e6, V);
fprintf('%-10s %10s %8s %8s %8s %8s\n', 'type', 'f0 re bi', '2f0', '3f0', '5f0', '7f0');
h1bi = dft1(X{2}, f0);
for k = 1:5
  h1 = dft1(X{k}, f0); r = zeros(1,4); nn = [2 3 5 7];
  for j = 1:4, r(j) = max(-99, 20*log10(dft1(X{k}, nn(j)*f0)/h1 + 1e-12)); end
  fprintf('%-10s %8.2f dB %7.1f %8.1f %8.1f %8.1f   dBc\n', names{k}, 20*log10(h1/h1bi), r);
end
fprintf('theory: bipolar 3f0 %.2f, 5f0 %.2f dBc; 3-level f0 %.2f dB; 5-level f0 %.2f dB\n', ...
  20*log10(1/3), 20*log10(1/5), 20*log10(sind(beta/2)), 20*log10(0.5*(cosd(12)+cosd(48))));

%% Part 2: slew-limited edges (constant-current charging of Cload)
fprintf('\nPart 2: slew limit I/C = %.2f V/ns, full -V to +V edge = %.1f ns\n', slew/1e9, 2*V*Cload/Ipk*1e9);
st = slew*dt; Y = cell(5,1);
for k = 1:5
  x = X{k}; y = zeros(size(x)); prev = 0;
  for i = 1:N
    d = x(i) - prev; prev = prev + max(-st, min(st, d)); y(i) = prev;
  end
  Y{k} = y;
  h1 = dft1(y, f0); r = zeros(1,3); nn = [3 5 7];
  for j = 1:3, r(j) = max(-99, 20*log10(dft1(y, nn(j)*f0)/h1 + 1e-12)); end
  fprintf('%-10s f0 %6.2f dB re ideal bipolar   3f0 %6.1f  5f0 %6.1f  7f0 %6.1f dBc\n', names{k}, 20*log10(h1/h1bi), r);
end
fprintf('trapezoid prediction for bipolar: f0 x sinc = %.2f dB, 3f0 rel. = %.1f dBc\n', ...
  20*log10(sin(pi*f0*2*V*Cload/Ipk)/(pi*f0*2*V*Cload/Ipk)), ...
  20*log10((1/3)*sin(3*pi*f0*2*V*Cload/Ipk)/(3*pi*f0*2*V*Cload/Ipk)/(sin(pi*f0*2*V*Cload/Ipk)/(pi*f0*2*V*Cload/Ipk))));

%% Part 3: duty error and the second harmonic (bipolar, ideal edges)
fprintf('\nPart 3: duty error, bipolar, 20 cycles\n');
kList = [1 2 4 7 11];                  % positive half-cycle longer by k samples
dl = zeros(size(kList)); sim2 = dl; thr2 = dl;
inb = (ph >= 0) & (ph < 20); fr = ph - floor(ph);
for j = 1:numel(kList)
  e  = (kList(j) - 0.5)/720;           % boundary placed between samples
  x  = V * lvBi(dutyMap(fr, e)) .* inb;
  npos = sum(x(inb) > 0)/20;           % realized samples per positive half cycle
  dl(j)   = (npos - 360)*dt;           % realized Delta (s)
  sim2(j) = 20*log10(dft1(x, 2*f0)/dft1(x, f0));
  thr2(j) = 20*log10(sin(pi*dl(j)*f0));
  fprintf('Delta = %5.2f ns   2f0 simulated %6.2f dBc   theory 20log10(sin(pi*Delta*f0)) %6.2f dBc\n', dl(j)*1e9, sim2(j), thr2(j));
end
fprintf('1 ns at 5 MHz (theory): %.2f dBc; a -45 dBc spec needs Delta < %.2f ns\n', ...
  20*log10(sin(pi*1e-9*f0)), asin(10^(-45/20))/(pi*f0)*1e9);

%% Figures
figure('Name','m23 pulser waveforms and spectra', 'Position', [50 50 900 950]);
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'PaperPosition', [0 0 10 11]);
nf = 2^16; fax = (0:nf/2-1)/(nf*dt)/1e6;
for k = 1:5
  y = Xs{k};
  ys = zeros(size(y)); prev = 0;
  for i = 1:N, d = y(i) - prev; prev = prev + max(-st, min(st, d)); ys(i) = prev; end
  subplot(5,2,2*k-1); idx = t < 1.2e-6;
  plot(t(idx)*1e6, y(idx), ':', t(idx)*1e6, ys(idx), '-'); ylim([-1.2 1.2]*V);
  ylabel('V'); title([names{k} ', 4 cycles']); if k == 5, xlabel('time (\mus)'); end
  S = abs(fft(ys, nf)); S = S(1:nf/2); S = 20*log10(S/max(S) + 1e-9);
  subplot(5,2,2*k); plot(fax, S); xlim([0 6*f0/1e6]); ylim([-70 3]); grid on;
  ylabel('dB'); title('spectrum (slew-limited)'); if k == 5, xlabel('frequency (MHz)'); end
end
figure('Name','m23 duty error');
dd = linspace(0.05e-9, 4e-9, 200);
plot(dd*1e9, 20*log10(sin(pi*dd*f0)), '-', dl*1e9, sim2, 'o'); grid on;
xlabel('duty error \Delta (ns)'); ylabel('2f_0 level (dBc)'); legend('sin(\pi\Delta f_0)', 'simulated', 'Location', 'southeast');
title('Second harmonic from a duty error, 5 MHz bipolar');
