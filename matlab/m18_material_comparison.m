%% m18_material_comparison.m
% Transducer materials side by side: bulk (d33, kt, eps, Z, temperature limit)
% and thin films (e31,f and the PMUT figure of merit e31,f^2 / (eps0 eps33)).
%
% What it shows: why single crystals win on coupling, composites on impedance,
%   PVDF on receive voltage (g33) and water match, and ScAlN/PZT films on FOM.
% What to change: add your own supplier values to the tables.
% Sources (see the Module 18 table): PZT-5H: Ma et al. IEEE TUFFC 2015 Table I
%   (d33 593, kt 0.51, epsS 1470, Z 36.0) and APC/Piezo Systems PSI-5H4E (k33 0.75,
%   epsT 3800, Tc 230 C); PMN-PT and PIN-PMN-PT: TRS Technologies X2B/X4B datasheet
%   ranges; composites: Smith-Auld model (m18_composite_effective_properties.m);
%   KNN: Saito et al. Nature 2004; PVDF: Measurement Specialties data (d33 -33,
%   g33 -330e-3, kt 0.14) and Precision Acoustics (rho 1780, c 2250, eps 8-10);
%   films: Roy, Lee & Lee, Microsyst. Nanoeng. 2023, Table I.
% Expected output: g33 (1e-3 V m/N) PZT-5H 18.5, PMN-PT 36.7, PIN-PMN-PT 37.3,
%   PVDF 356 (datasheet 330); film FOM (GPa) AlN 11.9, Sc0.15AlN 24.1,
%   sol-gel PZT 5.6-54.2, epitaxial PZT 74-111, KNN film 18-53.
% Runtime: < 1 s.
clear; close all; rng(1);
e0 = 8.854e-12;

%% Bulk materials: [low high] ranges; NaN = not given by the cited sources
names = {'PZT-5H', 'PMN-PT', 'PIN-PMN-PT', '1-3 PZT/epoxy', '1-3 crystal/epoxy', 'KNN (textured)', 'PVDF'};
d33  = [593 650; 1500 2400; 1200 2100; NaN NaN; NaN NaN; 300 416; 30 33];        % pC/N (PVDF magnitude)
kt   = [0.51 0.51; 0.55 0.59; 0.52 0.58; 0.67 0.69; 0.81 0.83; NaN NaN; 0.14 0.14];
epsT = [3800 3800; 5000 7000; 4000 6000; NaN NaN; NaN NaN; NaN NaN; 8 12];
Z    = [36.0 36.0; 36.4 39.2; 35.2 37.8; 11.0 21.9; 11.3 22.8; NaN NaN; 4.0 4.0];  % MRayl
Tmax = [230 230; 86 100; 105 130; NaN NaN; NaN NaN; NaN NaN; 70 80];             % C: Tc (PZT), Trt (crystals), max use (PVDF)
mid = @(x) mean(x, 2);
fprintf('%-18s %10s %6s %8s %8s %8s %12s\n', 'material', 'd33 pC/N', 'kt', 'epsT', 'Z MRayl', 'Tmax C', 'g33 mVm/N');
g33 = mid(d33)*1e-12./(mid(epsT)*e0)*1e3;
for i = 1:numel(names)
  fprintf('%-18s %10.0f %6.2f %8.0f %8.1f %8.0f %12.1f\n', names{i}, mid(d33(i,:)), mid(kt(i,:)), mid(epsT(i,:)), mid(Z(i,:)), mid(Tmax(i,:)), g33(i));
end
fprintf('PVDF g33 from the datasheet: 330 mV m/N (computed above from d33 and eps: %.0f)\n', g33(end));
fprintf('Intensity transmission into 1.5 MRayl at the mid impedance: ');
T = 4*mid(Z)*1.5./(mid(Z)+1.5).^2; fprintf('%s ', sprintf('%.2f', T(1))); fprintf('(PZT), %.2f (composite), %.2f (PVDF)\n', T(4), T(7));

%% Thin films for PMUTs (Roy et al. 2023, Table I)
fn  = {'AlN', 'Sc0.15AlN', 'PZT (sol-gel)', 'PZT (epitaxial)', 'KNN film'};
e31 = [1.05 1.05; 1.6 1.6; 8 12; 14 14; 8.5 14.4];      % |e31,f| C/m^2
epf = [10.5 10.5; 12 12; 300 1300; 200 300; 445 445];   % eps33,f
fprintf('\n%-16s %12s %10s %14s %14s\n', 'film', '|e31,f|', 'eps33', 'e31/eps0eps GV/m', 'e31^2/eps0eps GPa');
fomT = zeros(numel(fn),2);
for i = 1:numel(fn)
  rx = e31(i,:)./(e0*epf(i,[2 1]))/1e9;                 % range: low e31 with high eps, high e31 with low eps
  fom = e31(i,:).^2./(e0*epf(i,[2 1]))/1e9; fomT(i,:) = fom;
  fprintf('%-16s %5.2f-%5.2f %4.0f-%4.0f %6.1f-%6.1f %6.1f-%6.1f\n', fn{i}, e31(i,1), e31(i,2), epf(i,1), epf(i,2), rx(1), rx(2), fom(1), fom(2));
end

%% Figure
figure('Position', [80 80 1100 380]);
subplot(1,3,1); barh(mid(kt)); set(gca, 'YTick', 1:numel(names), 'YTickLabel', names); xlabel('k_t'); title('Thickness coupling'); grid on;
subplot(1,3,2); barh(mid(Z)); set(gca, 'YTick', 1:numel(names), 'YTickLabel', names); xlabel('Z (MRayl)'); title('Acoustic impedance'); grid on;
subplot(1,3,3); barh(mean(fomT,2)); set(gca, 'YTick', 1:numel(fn), 'YTickLabel', fn); xlabel('e_{31,f}^2/(\epsilon_0\epsilon_{33}) (GPa)'); title('Film FOM (range mid)'); grid on;
