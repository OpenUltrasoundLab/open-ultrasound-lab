%% m06_noise_floor_mp.m  --  Module 06: the noise floor of the SVD
% Purpose : Singular values of a pure-noise matrix fill the Marchenko-Pastur
%           bulk  sigma*(sqrt(m)-sqrt(n)) ... sigma*(sqrt(m)+sqrt(n)).
%           A rank-1 coherent component of strength theta separates from the
%           bulk only if theta > sigma*(m*n)^(1/4) (BBP transition).
% Change  : m, n (Casorati shape: pixels x frames), sig (noise std),
%           ratios (spike strength / threshold).
% Expect  : measured noise edges within ~2 % of theory for m = 1600, n = 80;
%           sigma_1 and singular-vector overlaps match the BBP formulas.
% Runtime : ~3 s in Octave.
clear; close all; rng(1);

%% 1. Pure complex Gaussian noise: bulk edges
sig = 1.0;                                  % noise std per entry, E|z|^2 = sig^2
shapes = [1600 80; 1600 400; 4096 128; 96 96];
fprintf('--- pure noise, complex Gaussian, sigma = %.1f\n', sig);
fprintf('   m      n   max(s) theory  min(s) theory  spread(dB)\n');
for i = 1:size(shapes,1)
    m = shapes(i,1); n = shapes(i,2);
    X = sig*(randn(m,n) + 1i*randn(m,n))/sqrt(2);
    s = svd(X);
    fprintf('%5d %5d  %6.2f %6.2f  %6.2f %6.2f   %5.1f\n', m, n, s(1), sig*(sqrt(m)+sqrt(n)), ...
        s(end), sig*(sqrt(m)-sqrt(n)), 20*log10(s(1)/s(end)));
end
fprintf('Gavish-Donoho optimal hard threshold (square, known sigma): %.3f*sqrt(n)*sigma (edge is 2*sqrt(n)*sigma)\n', 4/sqrt(3));

%% 2. Histogram of singular values against the Marchenko-Pastur density
m = 1600; n = 400; c = n/m;
X = sig*(randn(m,n) + 1i*randn(m,n))/sqrt(2);
s = svd(X);
a = (1-sqrt(c))^2; b = (1+sqrt(c))^2;      % eigenvalue support of X'X/(sig^2 m)
sg = linspace(sig*sqrt(m*a), sig*sqrt(m*b), 400);
x = sg.^2/(sig^2*m);
p = sqrt(max((b-x).*(x-a),0))./(2*pi*c*x) .* 2.*sg/(sig^2*m);   % density of s
figure('Position',[100 100 900 340]);
subplot(1,2,1);
[cnt, ctr] = hist(s, 30); bw = ctr(2)-ctr(1);
bar(ctr, cnt/(sum(cnt)*bw), 1, 'FaceColor', [0.75 0.8 0.9]); hold on;
plot(sg, p, 'r', 'LineWidth', 2);
xlabel('singular value'); ylabel('density');
title(sprintf('Pure noise %dx%d: MP density', m, n));
legend('svd(X)', 'Marchenko-Pastur', 'Location', 'north');

%% 3. Rank-1 coherent component in noise: BBP detectability threshold
m = 1600; n = 80; c = n/m;
thc = sig*(m*n)^(1/4);                     % detection threshold on theta
edge = sig*(sqrt(m)+sqrt(n));
fprintf('\n--- rank-1 spike in %dx%d noise: theta_c = %.2f, bulk edge = %.2f\n', m, n, thc, edge);
fprintf('theta/theta_c  s1(meas) s1(pred)  |<v1,v>|^2 meas/pred  |<u1,u>|^2 meas/pred\n');
ratios = [0.5 0.8 1.0 1.2 1.5 2 3 4];
ntr = 10; s1m = zeros(size(ratios));
for i = 1:numel(ratios)
    th = ratios(i)*thc; acc = zeros(1,3);
    for r = 1:ntr
        u = randn(m,1) + 1i*randn(m,1); u = u/norm(u);   % spatial pattern
        v = randn(n,1) + 1i*randn(n,1); v = v/norm(v);   % time course
        Y = th*(u*v') + sig*(randn(m,n) + 1i*randn(m,n))/sqrt(2);
        [U,S,V] = svd(Y, 'econ');
        acc = acc + [S(1,1), abs(V(:,1)'*v)^2, abs(U(:,1)'*u)^2]/ntr;
    end
    tp = th/(sig*sqrt(m));                 % normalized strength
    if ratios(i) > 1
        s1p = sqrt((th^2 + sig^2*m)*(th^2 + sig^2*n))/th;
        pv = 1 - c*(1+tp^2)/(tp^2*(tp^2+c)); pu = 1 - (c+tp^2)/(tp^2*(tp^2+1));
    else
        s1p = edge; pv = 0; pu = 0;
    end
    s1m(i) = acc(1);
    fprintf('   %4.1f        %6.2f   %6.2f      %5.2f / %5.2f          %5.2f / %5.2f\n', ...
        ratios(i), acc(1), s1p, acc(2), pv, acc(3), pu);
end
subplot(1,2,2);
tt = linspace(0.05, 4.2, 200)*thc;
s1th = edge*ones(size(tt)); k = tt > thc;
s1th(k) = sqrt((tt(k).^2 + sig^2*m).*(tt(k).^2 + sig^2*n))./tt(k);
plot(tt/thc, s1th, 'r', 'LineWidth', 2); hold on;
plot(ratios, s1m, 'ko', 'MarkerFaceColor', 'k');
plot([0 4.2], [edge edge], 'b--');
xlabel('\theta / \theta_c'); ylabel('\sigma_1');
title('Spike emerges only above \theta_c = \sigma (mn)^{1/4}');
legend('BBP prediction', 'measured', 'bulk edge', 'Location', 'northwest');
fprintf('Note: theta_c = %.1f sigma is far below the bulk edge %.1f sigma: a coherent\n', thc/sig, edge/sig);
fprintf('component is detectable even when sigma_1 barely clears the noise bulk.\n');
