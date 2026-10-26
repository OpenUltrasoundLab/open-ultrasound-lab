%% m06_svd_geometry.m  --  Module 06: SVD geometry, eig vs svd, power iteration, DORT
% Purpose : (1) the unit circle maps to an ellipse with semi-axes sigma_i*u_i;
%           (2) eigenvalues and singular values agree only for normal matrices;
%           (3) power iteration converges as |lambda2/lambda1|^k;
%           (4) DORT: the SVD of an array transfer matrix gives one singular
%               value per point scatterer and a focusing law for each.
% Change  : A (2x2), the gap ratio r, scatterer positions/reflectivities.
% Expect  : sigma = 1.618, 0.618 for the shear; power-iteration slope = log10(r);
%           two significant singular values for two point scatterers.
% Runtime : ~2 s in Octave.
clear; close all; rng(1);

%% 1. Circle -> ellipse
A = [1.2 0.9; -0.3 0.6];
[U,S,V] = svd(A);
th = linspace(0, 2*pi, 361); C = [cos(th); sin(th)]; E = A*C;
fprintf('A = [%.2f %.2f; %.2f %.2f]\n', A');
fprintf('sigma = %.4f, %.4f   |det A| = %.4f = sigma1*sigma2 = %.4f\n', S(1,1), S(2,2), abs(det(A)), S(1,1)*S(2,2));
fprintf('||A V - U S|| = %.1e,  ||A||_F^2 = %.4f = sum sigma^2 = %.4f\n', norm(A*V - U*S), norm(A,'fro')^2, sum(diag(S).^2));
fprintf('v1 at %.1f deg, u1 at %.1f deg (signs of (u_i, v_i) pairs are arbitrary)\n', ...
    atan2(V(2,1),V(1,1))*180/pi, atan2(U(2,1),U(1,1))*180/pi);
figure('Position',[100 100 900 420]);
subplot(1,2,1);
plot(C(1,:), C(2,:), 'Color', [0.6 0.6 0.6]); hold on; plot(E(1,:), E(2,:), 'b', 'LineWidth', 2);
for i = 1:2
    plot([0 V(1,i)], [0 V(2,i)], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
    plot([0 S(i,i)*U(1,i)], [0 S(i,i)*U(2,i)], 'r', 'LineWidth', 2);
end
axis equal; axis([-1.8 1.8 -1.8 1.8]); grid on; xlabel('x_1'); ylabel('x_2');
title('unit circle (grey) -> A x (blue); v_i (dashed) -> \sigma_i u_i (red)');

%% 2. Eigenvalues versus singular values
names = {'symmetric', 'non-symmetric', 'shear', 'rotation 35 deg x 1.1'};
t = 35*pi/180;
mats = {[1.5 0.5; 0.5 0.8], [1.2 0.9; 0.15 0.5], [1 1; 0 1], 1.1*[cos(t) -sin(t); sin(t) cos(t)]};
fprintf('\n%-24s %-26s %-16s %s\n', 'matrix', 'eig', 'svd', 'angle(v1,v2)');
for i = 1:numel(mats)
    M = mats{i}; [X,L] = eig(M); l = diag(L); s = svd(M);
    if isreal(l)
        ang = acos(min(1, abs(X(:,1)'*X(:,2))/(norm(X(:,1))*norm(X(:,2)))))*180/pi;
        fprintf('%-24s %8.4f %8.4f          %6.4f %6.4f   %5.1f deg\n', names{i}, l(1), l(2), s(1), s(2), ang);
    else
        fprintf('%-24s %6.3f%+6.3fj (pair)       %6.4f %6.4f   complex\n', names{i}, real(l(1)), abs(imag(l(1))), s(1), s(2));
    end
end

%% 3. Power iteration on an 8x8 symmetric matrix
n = 8; [Q,~] = qr(randn(n));
subplot(1,2,2);
rs = [0.5 0.7 0.95]; cols = {'b-', 'r--', 'k-.'}; hp = zeros(1, numel(rs));
for j = 1:numel(rs)
    r = rs(j);
    lam = [1, r, r*0.75.^(1:n-2)];
    A8 = Q*diag(lam)*Q';
    x = randn(n,1); K = 80; err = zeros(1,K+1);
    for k = 0:K
        x = x/norm(x); c1 = Q(:,1)'*x;
        err(k+1) = norm(x - c1*Q(:,1))/abs(c1);          % tan of angle to q1
        x = A8*x;
    end
    hp(j) = semilogy(0:K, err, cols{j}, 'LineWidth', 1.6); hold on;
    kk = 5:min(K, find(err > 1e-13, 1, 'last'));
    p = polyfit(kk, log10(err(kk+1)), 1);
    fprintf('power iteration r = %.2f: slope %.4f (theory %.4f), iterations to 1e-15 = %d\n', ...
        r, p(1), log10(r), ceil(-15/log10(r)));
end
ylim([1e-16 1e2]); grid on; xlabel('iteration k'); ylabel('tan \angle(x_k, q_1)');
title('Power iteration: error ~ r^k'); legend(hp, 'r = 0.5', 'r = 0.7', 'r = 0.95', 'Location', 'southwest');

%% 4. DORT: SVD of the inter-element transfer matrix
c0 = 1540; f0 = 5e6; lam0 = c0/f0; kw = 2*pi/lam0;   % m, rad/m
Ne = 64; pitch = 0.3e-3; xe = ((0:Ne-1) - (Ne-1)/2)*pitch;  % element x (m)
sc = [-2e-3 30e-3 1.0; 3e-3 32e-3 0.5];               % [x z reflectivity] of two scatterers
G = zeros(Ne, size(sc,1));
for q = 1:size(sc,1)
    d = sqrt((xe - sc(q,1)).^2 + sc(q,2)^2);
    G(:,q) = exp(-1i*kw*d(:))./sqrt(d(:));             % 2-D Green's function (far-field form)
end
K = G*diag(sc(:,3))*G.';                               % K = K.' by reciprocity
K = K + 1e-3*norm(K)*(randn(Ne) + 1i*randn(Ne))/sqrt(2*Ne);
[Uk,Sk,Vk] = svd(K); sk = diag(Sk);
fprintf('\nDORT, 64 elements, two scatterers (refl. 1.0, 0.5): sigma_1..4 / sigma_1 = %.3f %.3f %.4f %.4f\n', sk(1:4)/sk(1));
fprintf('lateral resolution lambda*z/D = %.2f mm; lateral separation of scatterers = %.1f mm\n', lam0*31e-3/(Ne*pitch)*1e3, (sc(2,1)-sc(1,1))*1e3);
figure('Position',[100 100 900 340]);
subplot(1,2,1); semilogy(sk/sk(1), 'ko-'); grid on; xlabel('index i'); ylabel('\sigma_i / \sigma_1');
title('DORT: one singular value per scatterer');
xg = linspace(-8e-3, 8e-3, 401); zf = 31e-3;          % back-propagate onto a line at 31 mm
subplot(1,2,2); hq = zeros(1,2);
for q = 1:2
    w = conj(Uk(:,q));                                  % phase-conjugate (time-reversed) law
    P = zeros(size(xg));
    for ix = 1:numel(xg)
        d = sqrt((xe - xg(ix)).^2 + zf^2);
        P(ix) = abs(sum(w(:).*exp(-1i*kw*d(:))./sqrt(d(:))));
    end
    hq(q) = plot(xg*1e3, 20*log10(P/max(P)), 'LineWidth', 1.6); hold on;
    [~, im] = max(P); fprintf('singular vector %d focuses at x = %.2f mm\n', q, xg(im)*1e3);
end
hq(3) = plot(sc(:,1)*1e3, [0 0], 'kv', 'MarkerFaceColor', 'k');
ylim([-40 3]); grid on; xlabel('lateral position (mm)'); ylabel('field (dB)');
title('Back-propagated u_1, u_2 at z = 31 mm'); legend(hq, 'conj(u_1)', 'conj(u_2)', 'scatterers', 'Location', 'northwest');
