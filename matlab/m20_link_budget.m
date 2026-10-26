%% m20_link_budget.m
% Link budget from transmit voltage to beamformed SNR (Module 20.8, Demo 20.7).
%   Voc   = Vtx * 10^((Gloop + R - 2*alpha*f0*z)/20)     open-circuit echo (peak)
%   H, Zs from the BVD element + cable + LNA model at f0 (as m20_noise_budget.m)
%   noise at the LNA input: (4kT Re Zs + en^2 + in^2 |Zs|^2) * B, B = fbw*f0
%   TGC: VGA gain puts the echo at the target dBFS, clamped to [vmin, vmax]
%   ADC noise in band: full-scale-sine rms * 10^(-SNRadc/20) * sqrt(B/(fs/2))
%   SNRout = SNRch + 10log10(Gbf), Gbf = N (point) or (2N^2+1)/(3N) (speckle)
% What to change: any parameter in the first cell (drive, f0, attenuation,
%   cable, LNA, VGA range, bits, channels, target).
% Expected output (defaults): Voc -96.02 dBV peak, |H| -12.00 dB, input noise
%   2.06 uV rms, SNRch 2.62 dB, VGA +30 dB, -58.0 dBFS, ADC penalty 0.07 dB,
%   SNRout 21.93 dB, penetration 10.4 cm; without the cable 24.83 dB.
clear; close all; rng(1);

%% Parameters (the Demo 20.7 defaults)
P.Vtx  = 50;      % drive voltage, peak (V)
P.f0   = 5e6;     % centre frequency (Hz)
P.fbw  = 0.7;     % noise bandwidth as a fraction of f0
P.al   = 0.5;     % attenuation (dB/cm/MHz)
P.z    = 6;       % depth (cm)
P.Gl   = -40;     % element pulse-echo loop gain, perfect reflector at focus (dB), illustrative
P.R    = -60;     % target reflectivity re perfect reflector (dB): speckle
P.speck= 1;       % 1 = speckle coherence gain (2N^2+1)/(3N), 0 = point target (N)
P.Cc   = 150e-12; % cable capacitance (F)
P.en   = 0.8e-9;  % LNA voltage noise (V/rtHz)
P.in   = 1.5e-12; % LNA current noise (A/rtHz)
P.Glna = 20;      % LNA gain (dB)
P.vmin = -10; P.vmax = 30;   % VGA gain range (dB)
P.tgt  = -40;     % TGC target for the echo (dBFS)
P.bits = 12; P.Limp = 3;     % ADC bits, implementation loss (dB)
P.fs   = 40e6;    % sample rate (Hz)
P.Vfs  = 2;       % ADC full scale (Vpp)
P.N    = 128;     % receive channels
kB = 1.380649e-23; T0 = 290;
% element model (BVD) and input capacitance
C0 = 30e-12; k2 = 0.4; Q = 2; Cin = 15e-12;

%% Model (anonymous functions of the parameter struct p)
w0f  = @(p) 2*pi*p.f0;
Cmf  = C0*k2/(1-k2);
Zef  = @(p) 1 ./ (1./(1/(w0f(p)*Cmf*Q)) + 1i*w0f(p)*C0);   % at resonance the motional branch is Rm
Zlf  = @(p) 1 ./ (1i*w0f(p)*(p.Cc + Cin));
Hf   = @(p) abs(Zlf(p) ./ (Zef(p) + Zlf(p)));
Zsf  = @(p) Zef(p).*Zlf(p) ./ (Zef(p) + Zlf(p));
Sn   = @(p) 4*kB*T0*real(Zsf(p)) + p.en^2 + p.in^2*abs(Zsf(p))^2;     % V^2/Hz at LNA input
Bf   = @(p) p.fbw*p.f0;
Vocf = @(p) p.Vtx * 10^((p.Gl + p.R - 2*p.al*p.f0/1e6*p.z)/20);
Vinf = @(p) Vocf(p)*Hf(p);
vnf  = @(p) sqrt(Sn(p)*Bf(p));
gf   = @(p) min(p.vmax, max(p.vmin, p.tgt - 20*log10(Vinf(p)*10^(p.Glna/20)/(p.Vfs/2))));
Gtf  = @(p) p.Glna + gf(p);
vqf  = @(p) (p.Vfs/2/sqrt(2)) * 10^(-(6.02*p.bits + 1.76 - p.Limp)/20) * sqrt(Bf(p)/(p.fs/2));
snrch= @(p) 20*log10( (Vinf(p)*10^(Gtf(p)/20)/sqrt(2)) / sqrt((vnf(p)*10^(Gtf(p)/20))^2 + vqf(p)^2) );
gbf  = @(p) p.speck*10*log10((2*p.N^2+1)/(3*p.N)) + (1-p.speck)*10*log10(p.N);
snro = @(p) snrch(p) + gbf(p);

%% Stage table at the default depth
p = P; dBv = @(v) 20*log10(v);
Gt = Gtf(p); Va = Vinf(p)*10^(Gt/20); vna = vnf(p)*10^(Gt/20);
fprintf('%-34s %9s %12s %12s %9s\n', 'stage', 'gain dB', 'signal dBV', 'noise dBV', 'SNR dB');
fprintf('%-34s %9s %12.2f\n', sprintf('drive %g V peak (rms)', p.Vtx), '', dBv(p.Vtx/sqrt(2)));
fprintf('%-34s %9.2f\n', 'loop gain', p.Gl);
fprintf('%-34s %9.2f\n', 'target', p.R);
fprintf('%-34s %9.2f\n', 'two-way attenuation', -2*p.al*p.f0/1e6*p.z);
elN = sqrt(4*kB*T0*real(Zef(p))*Bf(p));
fprintf('%-34s %9s %12.2f %12.2f %9.2f\n', 'open-circuit echo (rms)', '', dBv(Vocf(p)/sqrt(2)), dBv(elN), dBv(Vocf(p)/sqrt(2)/elN));
fprintf('%-34s %9.2f\n', 'cable / interface |H|', dBv(Hf(p)));
fprintf('%-34s %9s %12.2f %12.2f %9.2f\n', 'LNA input', '', dBv(Vinf(p)/sqrt(2)), dBv(vnf(p)), dBv(Vinf(p)/sqrt(2)/vnf(p)));
fprintf('%-34s %9.2f\n', 'LNA gain', p.Glna);
fprintf('%-34s %9.2f\n', 'VGA (TGC) gain', gf(p));
fprintf('%-34s %9s %12.2f %12.2f %9.2f\n', 'ADC input', '', dBv(Va/sqrt(2)), dBv(sqrt(vna^2 + vqf(p)^2)), snrch(p));
fprintf('%-34s %9.2f\n', 'beamforming gain', gbf(p));
fprintf('%-34s %9s %12s %12s %9.2f\n', 'output SNR', '', '', '', snro(p));
fprintf('\nVoc = %.2f dBV peak (%.2f uV), Vin = %.2f uV peak, input noise %.3f uV rms (%.3f nV/rtHz)\n', ...
  dBv(Vocf(p)), Vocf(p)*1e6, Vinf(p)*1e6, vnf(p)*1e6, sqrt(Sn(p))*1e9);
fprintf('echo at ADC %.1f dBFS, analog noise at ADC %.3f mV rms, ADC noise in band %.1f uV rms\n', dBv(Va/(p.Vfs/2)), vna*1e3, vqf(p)*1e6);
q = p; q.Limp = -300; pen = snrch(q) - snrch(p);
fprintf('ADC penalty %.3f dB; LSB %.1f uV; analog noise = %.2f LSB rms\n', pen, p.Vfs/2^p.bits*1e6, vna/(p.Vfs/2^p.bits));

%% SNR versus depth and penetration
z = 0.1:0.01:25; s = zeros(size(z)); sa = s; sq = s;
for j = 1:numel(z)
  q = P; q.z = z(j); s(j) = snro(q);
  q.Limp = -300; sa(j) = snro(q);                      % ideal ADC
end
iz = find(s(1:end-1) > 0 & s(2:end) <= 0, 1);
zpen = z(iz) + 0.01*s(iz)/(s(iz) - s(iz+1));
fprintf('\nnoise-limited penetration (SNRout = 0 dB): %.2f cm\n', zpen);
q = P; q.Cc = 0;       fprintf('no cable: SNRout %.2f dB (+%.2f dB)\n', snro(q), snro(q) - snro(P));
q = P; q.Vtx = 100;    fprintf('100 V drive: SNRout %.2f dB\n', snro(q));
q = P; q.R = -20; q.speck = 0; fprintf('strong interface (point) at 6 cm: SNRout %.2f dB, echo %.1f dBFS\n', snro(q), dBv(Vinf(q)*10^(Gtf(q)/20)/(q.Vfs/2)));
q = P; q.vmax = 0;     fprintf('VGA max 0 dB: SNRout %.2f dB (ADC-limited)\n', snro(q));

%% Figures
figure('Name','m23 link budget');
st = [dBv(P.Vtx/sqrt(2)), P.Gl, P.R, -2*P.al*P.f0/1e6*P.z, dBv(Hf(P)), P.Glna, gf(P)];
lev = cumsum(st);
subplot(2,1,1);
bar(1:7, lev); hold on; plot([4 5 6 7], dBv([elN vnf(P) vnf(P)*10^(P.Glna/20) vna]), 'o-'); hold off; grid on;
set(gca, 'XTick', 1:7, 'XTickLabel', {'drive','loop','target','atten','cable','LNA','VGA'});
ylabel('dBV rms'); title('cumulative signal level (bars) and noise floor (line)');
subplot(2,1,2);
plot(z, s, '-', z, sa, '--', [P.z P.z], [-20 80], ':'); grid on; ylim([-20 80]);
xlabel('depth (cm)'); ylabel('SNR_{out} (dB)'); legend('full chain', 'ideal ADC', 'selected depth'); title('beamformed SNR vs depth, speckle');
