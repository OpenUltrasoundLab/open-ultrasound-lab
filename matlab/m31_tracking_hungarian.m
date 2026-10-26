%% m31_tracking_hungarian.m
% Purpose : frame-to-frame linking of ULM localizations. Compares nearest
%           neighbour (NN), greedy global assignment, the optimal (Hungarian)
%           assignment and Hungarian with a constant-velocity prediction (the
%           predict step of a Kalman filter), against the Poisson prediction
%           for NN: P(correct) = exp(-pi * n * d^2).
% Model   : N bubbles on a unit torus (no edges), Poisson positions, density
%           n = N. Each bubble moves by a displacement of fixed length d in a
%           random direction; between frames its velocity turns by a random
%           vector of rms sv*d; localizations carry Gaussian errors sl*d per
%           axis. No births or deaths (every bubble has a partner).
%           Dimensionless crowding: rho = n * d^2 (bubbles per displacement^2).
% Hungarian: e-maxx / Jonker-style shortest augmenting path with potentials,
%           O(N^3), written out below (vectorized over columns), checked
%           against brute force on small problems.
% Change  : N, rhoList, nTrial, sl, sv.
% Expect  : NN follows exp(-pi rho) (0.729 measured vs 0.730 at rho = 0.1);
%           greedy global assignment is WORSE than NN (0.627): it locks in
%           short wrong pairs early. Hungarian is best without prediction
%           (0.856 at 0.1, 0.716 at 0.2). Prediction shrinks the residual to
%           0.40 d and keeps 0.96 correct at rho = 0.2 and 0.80 at rho = 1.
%           At 1 kHz, rho in ULM is 1e-5 to 2e-2 (table at the end), so in
%           practice links fail from missed or overlapping detections, not from
%           assignment. Runtime about 10 s (Octave).
rng(1);

N = 60; rhoList = [0.01 0.02 0.05 0.1 0.2 0.5 1 2]; nTrial = 40;
sl = 0.1; sv = 0.2;                                       % noise and velocity change, in units of d
names = {'NN', 'greedy', 'Hungarian', 'Hung.+pred'}; nM = 4;

%% Check the Hungarian code against brute force (6 x 6)
pm = perms(1:6); worst = 0;
for trial = 1:200
    C = rand(6);
    % ---- Hungarian (rows -> columns), minimal total cost ----
    n = 6; m = 6; u = zeros(n,1); v = zeros(1,m+1); p = zeros(1,m+1); way = zeros(1,m+1);
    for i = 1:n
        p(1) = i; j0 = 1; minv = inf(1,m+1); used = false(1,m+1);
        while true
            used(j0) = true; i0 = p(j0); cols = find(~used); cols = cols(cols > 1);
            cur = C(i0, cols-1) - u(i0) - v(cols); upd = cur < minv(cols);
            minv(cols(upd)) = cur(upd); way(cols(upd)) = j0;
            [delta, k] = min(minv(cols)); j1 = cols(k);
            uj = find(used); u(p(uj)) = u(p(uj)) + delta; v(uj) = v(uj) - delta; minv(~used) = minv(~used) - delta;
            j0 = j1; if p(j0) == 0, break; end
        end
        while true, j1 = way(j0); p(j0) = p(j1); j0 = j1; if j0 == 1, break; end, end
    end
    asg = zeros(1,n); asg(p(2:end)) = 1:m;
    costH = sum(C(sub2ind([n n], 1:n, asg)));
    costB = min(sum(C(sub2ind([n n], repmat(1:n, size(pm,1), 1), pm)), 2));
    worst = max(worst, abs(costH - costB));
end
fprintf('Hungarian versus brute force on 200 random 6 x 6 problems: max cost difference = %.2e\n', worst);

%% Linking experiment
res = zeros(numel(rhoList), nM); tic;
for ir = 1:numel(rhoList)
    d = sqrt(rhoList(ir)/N);                               % displacement for this crowding
    ok = zeros(1, nM); tot = 0;
    for trial = 1:nTrial
        x0 = rand(N,2);                                    % positions at frame t-1
        th = 2*pi*rand(N,1); v1 = d*[cos(th) sin(th)];     % velocity t-1 -> t
        v2 = v1 + sv*d/sqrt(2)*randn(N,2); v2 = d*v2./sqrt(sum(v2.^2,2));   % turned velocity t -> t+1, length d
        x1 = mod(x0 + v1, 1); x2 = mod(x1 + v2, 1);
        y0 = mod(x0 + sl*d*randn(N,2), 1); y1 = mod(x1 + sl*d*randn(N,2), 1); y2 = mod(x2 + sl*d*randn(N,2), 1);
        perm = randperm(N); y2 = y2(perm,:);               % unknown order in frame t+1; true partner of i is find(perm == i)
        truthIdx = zeros(N,1); truthIdx(perm) = 1:N;
        wrapd = @(a) a - round(a);                         % torus difference
        DX = wrapd(y1(:,1) - y2(:,1)'); DZ = wrapd(y1(:,2) - y2(:,2)'); D2 = DX.^2 + DZ.^2;
        vp = [wrapd(y1(:,1) - y0(:,1)) wrapd(y1(:,2) - y0(:,2))];   % velocity from the previous (correct) link
        PX = wrapd(y1(:,1) + vp(:,1) - y2(:,1)'); PZ = wrapd(y1(:,2) + vp(:,2) - y2(:,2)'); P2 = PX.^2 + PZ.^2;
        % NN: each bubble takes its nearest candidate (no exclusivity)
        [~, nnIdx] = min(D2, [], 2); ok(1) = ok(1) + sum(nnIdx == truthIdx);
        % greedy global: sort all pairs, accept if both are free
        [~, ord] = sort(D2(:)); ua = false(N,1); ub = false(N,1); gIdx = zeros(N,1); cnt = 0;
        for q = 1:numel(ord)
            [a, b] = ind2sub([N N], ord(q));
            if ~ua(a) && ~ub(b), ua(a) = true; ub(b) = true; gIdx(a) = b; cnt = cnt + 1; if cnt == N, break; end, end
        end
        ok(2) = ok(2) + sum(gIdx == truthIdx);
        % Hungarian on squared distances, then on squared prediction residuals
        for cm = 1:2
            if cm == 1, C = D2; else, C = P2; end
            n = N; m = N; u = zeros(n,1); v = zeros(1,m+1); p = zeros(1,m+1); way = zeros(1,m+1);
            for i = 1:n
                p(1) = i; j0 = 1; minv = inf(1,m+1); used = false(1,m+1);
                while true
                    used(j0) = true; i0 = p(j0); cols = find(~used); cols = cols(cols > 1);
                    cur = C(i0, cols-1) - u(i0) - v(cols); upd = cur < minv(cols);
                    minv(cols(upd)) = cur(upd); way(cols(upd)) = j0;
                    [delta, k] = min(minv(cols)); j1 = cols(k);
                    uj = find(used); u(p(uj)) = u(p(uj)) + delta; v(uj) = v(uj) - delta; minv(~used) = minv(~used) - delta;
                    j0 = j1; if p(j0) == 0, break; end
                end
                while true, j1 = way(j0); p(j0) = p(j1); j0 = j1; if j0 == 1, break; end, end
            end
            asg = zeros(N,1); asg(p(2:end)) = 1:m;
            ok(2 + cm) = ok(2 + cm) + sum(asg == truthIdx);
        end
        tot = tot + N;
    end
    res(ir,:) = ok/tot;
end
deff = sqrt(sv^2 + 12*sl^2);   % rms 2D prediction residual / d: velocity change plus noise of y(t+1) - 2y(t) + y(t-1) (1+4+1 per axis)
fprintf('\nFraction of correct links, N = %d per frame, %d trials, noise %.2f d, velocity change %.2f d rms\n', N, nTrial, sl, sv);
fprintf('   rho=n d^2  exp(-pi rho)  %s   1/(1+pi rho deff^2)\n', sprintf('%11s', names{:}));
for ir = 1:numel(rhoList)
    fprintf('   %6.2f     %8.3f    %s     %8.3f\n', rhoList(ir), exp(-pi*rhoList(ir)), sprintf('%11.3f', res(ir,:)), 1/(1 + pi*rhoList(ir)*deff^2));
end
fprintf('  (prediction residual: Gaussian, rms deff = sqrt(sv^2 + 12 sl^2) d = %.2f d, so NN with prediction gives 1/(1 + pi rho deff^2))\n', deff);
fprintf('Runtime %.1f s\n', toc);

%% What rho means for ULM numbers
FR = 1000; vb = [2 10 30];                                 % bubble speed, mm/s
nmm = [1 5 20];                                            % localizations per mm^2 per frame
fprintf('\nrho = n (v/FR)^2 at %d frames/s:\n   n (/mm^2)  v = %s mm/s\n', FR, sprintf('%6.0f', vb));
for k = 1:numel(nmm), fprintf('   %6.0f    %s\n', nmm(k), sprintf('%9.2e', nmm(k)*(vb/FR).^2)); end

%% Figure
figure('Position', [50 50 800 380]);
semilogx(rhoList, exp(-pi*rhoList), 'k-', 'LineWidth', 2); hold on;
mk = {'o-', 's-', 'd-', '^-'}; for m = 1:nM, semilogx(rhoList, res(:,m), mk{m}); end
semilogx(rhoList, 1./(1 + pi*rhoList*deff^2), 'k--');
xlabel('\rho = n d^2 (bubbles per displacement squared)'); ylabel('fraction of correct links');
legend(['e^{-\pi\rho} (NN theory)', names, '1/(1+\pi\rho d_{eff}^2)'], 'Location', 'southwest'); grid on; ylim([0 1.02]);
title('Linking accuracy versus crowding');
