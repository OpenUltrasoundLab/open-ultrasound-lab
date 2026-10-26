%% m30_lamb_a0.m
% Purpose : A0 (flexural) Lamb mode of an incompressible elastic plate in
%           vacuum, from the Rayleigh-Lamb antisymmetric equation with
%           c_L -> infinity:
%             (k^2 + q'^2)^2 tanh(k h/2) - 4 k^3 q' tanh(q' h/2) = 0,
%             k = w/c, q' = w sqrt(1/c^2 - 1/c_t^2), c < c_t.
%           Compared with the thin-plate limit sqrt(w h c_t / sqrt(3)), the
%           immersed-plate empirical fit sqrt(w h c_t / (2 sqrt(3))) of
%           Couade et al. 2010, and the incompressible Rayleigh speed.
% Change  : h (thickness), ct (bulk shear speed), frequency list.
% Expect  : Rayleigh speed 0.9553 c_t. h = 1 mm, c_t = 3 m/s: A0 = 1.365 m/s
%           at 200 Hz (apparent E 21 % of true), 2.675 m/s at 2 kHz.
rng(1);
Fn = @(c, f, h, ct) ((2*pi*f./c).^2 + (2*pi*f).^2.*max(1./c.^2 - 1/ct^2, 0)).^2.*tanh(pi*f*h./c) ...
     - 4*(2*pi*f./c).^3.*(2*pi*f.*sqrt(max(1./c.^2 - 1/ct^2, 0))).*tanh(pi*f*h*sqrt(max(1./c.^2 - 1/ct^2, 0)));
r = roots([1 -8 24 -16]); r = r(abs(imag(r)) < 1e-9); fprintf('Incompressible Rayleigh speed: %.4f c_t\n', sqrt(r));
cases = [1e-3 3; 1e-3 5; 0.55e-3 2]; flist = [50 100 200 500 1000 2000];
for ic = 1:size(cases,1)
  h = cases(ic,1); ct = cases(ic,2); fprintf('h = %.2f mm, c_t = %.1f m/s\n', h*1e3, ct);
  for f = flist
    cg = ct*linspace(1e-3, 0.999999, 600); v = Fn(cg, f, h, ct)./(2*pi*f./cg).^4;
    s = find(sign(v(1:end-1)) ~= sign(v(2:end)), 1, 'last'); a = cg(s); b = cg(s+1);
    g = @(c) Fn(c, f, h, ct)./(2*pi*f/c)^4; fa = g(a);
    for it = 1:60, mm = (a + b)/2; fm = g(mm); if sign(fm) == sign(fa), a = mm; fa = fm; else, b = mm; end, end
    cA = (a + b)/2; cK = sqrt(2*pi*f*h*ct/sqrt(3)); cF = sqrt(2*pi*f*h*ct/(2*sqrt(3)));
    fprintf('  %5d Hz: A0 = %.3f m/s, thin plate = %.3f, immersed fit = %.3f, (c/c_t)^2 = %.3f\n', f, cA, cK, cF, (cA/ct)^2);
  end
end
h = 1e-3; ct = 3; ff = 25:25:2000; cA = zeros(size(ff));
for i = 1:numel(ff)
  f = ff(i); cg = ct*linspace(1e-3, 0.999999, 600); v = Fn(cg, f, h, ct)./(2*pi*f./cg).^4; s = find(sign(v(1:end-1)) ~= sign(v(2:end)), 1, 'last');
  a = cg(s); b = cg(s+1); g = @(c) Fn(c, f, h, ct)./(2*pi*f/c)^4; fa = g(a);
  for it = 1:60, mm = (a + b)/2; fm = g(mm); if sign(fm) == sign(fa), a = mm; fa = fm; else, b = mm; end, end; cA(i) = (a + b)/2;
end
figure(1); plot(ff, cA, 'b', ff, sqrt(2*pi*ff*h*ct/sqrt(3)), 'r--', ff, sqrt(2*pi*ff*h*ct/(2*sqrt(3))), 'k:', ff, 0*ff + ct, 'k-', ff, 0*ff + 0.9553*ct, 'k--');
ylim([0 1.15*ct]); xlabel('f (Hz)'); ylabel('phase velocity (m/s)'); legend('A0 exact', 'thin plate', 'immersed (Couade)', 'c_t', 'Rayleigh', 'location', 'southeast'); title('1 mm plate, c_t = 3 m/s');
