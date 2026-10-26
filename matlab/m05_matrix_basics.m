%% m05_matrix_basics.m  Module 05: the matrix workbench (Demo 05.2) in numbers
% Purpose : print every operation of the Module 05 workbench on the same small
%           matrices as the page, so each number on the page can be checked here:
%           A*x two ways, A*B vs B*A, .' vs ' on IQ data, inv/det of 2x2 and 3x3,
%           A\b (square and tall), rank of a rank-deficient matrix, projectors,
%           the DFT matrix F4, convolution as a Toeplitz matrix, kron for
%           per-element operators on stacked channel data, and a tiny sparse DAS matrix.
% Change  : edit any matrix below; e.g. make row 3 of A = row 1 + row 2 and rerun.
% Expect  : A*x = [4 3 -2]; det(A) = 24; A*B = [2 1 3; 1 5 -5; 4 5 -1] ~= B*A;
%           z.'*z = 0 but z'*z = 2; A\b = [1 -1 0.5]; tall fit slope -0.98 dB/MHz;
%           rank(M) = 2 with null vector [1 1 -1]; F4'*F4 = 4*I; T*x = conv(h,x);
%           DAS: nnz(B) = 10, B*d = [2 0.75]; hilb(10): A\b residual ~1e-16 vs
%           inv(A)*b ~4e-5 at equal forward error ~5e-4 (backslash is backward stable).
% Runtime : < 1 s.  Base MATLAB only; runs in MATLAB R2019b+ and GNU Octave 8.
rng(1);
pm = @(name, M) fprintf(['%s =\n' repmat([repmat('%9.4g', 1, size(M,2)) '\n'], 1, size(M,1))], name, M.');   % print a real matrix row by row

%% 1. A*x: row picture (one dot product per output) and column picture (mix of columns)
A = [ 2 1 0;
     -1 3 2;
      0 1 4];                  % 3x3 real, not symmetric
x = [1; 2; -1];                % 3x1 column vector
for i = 1:3
    fprintf('row %d: [%g %g %g] . [%g %g %g] = %g\n', i, A(i,:), x, A(i,:)*x);
end
y = zeros(3,1);
for j = 1:3
    y = y + x(j)*A(:,j);       % add x_j times column j
    fprintf('after column %d: [%g %g %g]\n', j, y);
end
fprintf('A*x = [%g %g %g]\n\n', A*x);

%% 2. A*B versus B*A: the order matters
B = [1 0 2;
     0 1 -1;
     1 1 0];
pm('A*B', A*B);  pm('B*A', B*A);  pm('A*B - B*A', A*B - B*A);
fprintf('(A*B)(2,3) = %g*%g + %g*%g + %g*%g = %g\n\n', A(2,1), B(1,3), A(2,2), B(2,3), A(2,3), B(3,3), A(2,:)*B(:,3));

%% 3. Transpose .' versus conjugate transpose ' on IQ data
Z = [1+2i, 0.5-1i, -1i;
     2,   -1+1i,   0.3+0.4i];  % K = 2 snapshots (rows) x N = 3 elements (columns), complex IQ
disp('Z.'' (transpose, no conjugate) ='); disp(Z.');
disp('Z''  (conjugate transpose)    ='); disp(Z');
z = [1; 1i];
fprintf('z = [1; 1i]:  z.''*z = %g   z''*z = %g\n', z.'*z, z'*z);
R = Z'*Z;                      % 3x3 sample covariance across elements (times K): Hermitian PSD
disp('R = Z''*Z ='); disp(R);
fprintf('diag(R) = [%.2f %.2f %.2f] (element energies, real)\n', real(diag(R)));
ev = sort(real(eig(R)), 'descend');
fprintf('eig(R) = [%.4f %.4f %.4f], rank(R) = %d (K = 2 snapshots)\n\n', ev, rank(R));

%% 4. Inverse and determinant, 2x2 and 3x3
A2 = [4 7; 2 6];
fprintf('det(A2) = %g, inv(A2) = (1/det)[d -b; -c a]:\n', det(A2)); disp(inv(A2));
fprintf('det(A) = %g\n', det(A));
pm('24*inv(A) (the adjugate)', round(24*inv(A)));
pm('inv(A)', inv(A));
fprintf('max|A*inv(A) - I| = %.1e\n\n', max(max(abs(A*inv(A) - eye(3)))));

%% 5. Solve A x = b: square (LU) and tall (least squares by QR)
b = [1; -3; 1];
fprintf('A\\b = [%g %g %g]\n', A\b);
f  = [2; 4; 6; 8];                     % MHz
yd = [-1.0; -3.2; -4.8; -7.0];         % log-spectral ratio, dB
Af = [ones(4,1) f];                    % 4x2: intercept and slope columns
cf = Af \ yd;                          % least squares (QR inside)
r  = yd - Af*cf;
fprintf('tall 4x2 fit: c = Af\\y = [%.4f %.4f]  (intercept dB, slope dB/MHz)\n', cf);
fprintf('residual r = [%.3f %.3f %.3f %.3f], Af''*r = [%.1e %.1e]\n\n', r, Af'*r);

%% 6. Rank of a rank-deficient matrix
M = [1 2 3;
     2 0 2;
     0 1 1;
     1 1 2];                    % column 3 = column 1 + column 2
s = svd(M);
fprintf('svd(M) = [%.4f %.4f %.2e], rank(M) = %d\n', s, rank(M));
nv = [1; 1; -1];
fprintf('M*[1;1;-1] = [%g %g %g %g] (null vector)\n\n', M*nv);

%% 7. Projectors: P^2 = P
N = 4; one = ones(N,1);
P = eye(N) - one*one'/N;        % removes the slow-time mean (simplest wall filter)
pm('P = I - 11''/4', P);
xs = [3.0; 3.2; 2.9; 3.1];      % slow-time samples at one pixel: tissue (mean) + small flow
fprintf('P*x = [%.2f %.2f %.2f %.2f] (mean %.2f removed)\n', P*xs, mean(xs));
fprintf('max|P^2 - P| = %.1e, eig(P) = [%s], trace = %g\n', max(max(abs(P*P - P))), num2str(sort(eig(P)).', '%.0f '), trace(P));
H = Af * ((Af'*Af) \ Af');      % hat matrix of the line fit in section 5 (fine here: kappa small)
pm('H = Af*(Af''*Af)^-1*Af''', H);
fprintf('H*y = [%.3f %.3f %.3f %.3f] = Af*c, max|H^2 - H| = %.1e, trace(H) = %g\n\n', H*yd, max(max(abs(H*H - H))), trace(H));

%% 8. DFT matrix F4: F'*F = 4 I
n = (0:3);
F = exp(-2i*pi*n.'*n/4);        % F(k+1,n+1) = exp(-j 2 pi k n / 4): entries 1, -j, -1, j
disp('F4 ='); disp(round(F));   % round removes 1e-16 residue in printing
disp('F4''*F4 ='); disp(round(real(F'*F)*1e12)/1e12);
xf = [1; 2; 0; -1];
Fx = F*xf; Fx(abs(Fx) < 1e-12) = 0; Fx = round(Fx*1e12)/1e12;   % clean 1e-16 rounding residue for printing
fprintf('F*x    = %s\nfft(x) = %s\n', mat2str(Fx.', 4), mat2str(fft(xf).', 4));
fprintf('||F*x||^2 = %g = 4*||x||^2 = %g\n\n', norm(F*xf)^2, 4*norm(xf)^2);

%% 9. Convolution as a Toeplitz matrix
h  = [2; 1; -1];                % 3-tap pulse
xr = [1; 0; 0; 2; 0];           % reflectivity: two scatterers
T  = toeplitz([h; zeros(4,1)], [h(1) zeros(1,4)]);   % (5+3-1) x 5 = 7x5
pm('T = toeplitz([h; 0 0 0 0], [h(1) 0 0 0 0])', T);
fprintf('T*x       = [%s]\nconv(h,x) = [%s]\n\n', num2str((T*xr).', '%g '), num2str(conv(h, xr).', '%g '));

%% 10. kron: a per-element operator on stacked channel data d = X(:)
X  = [1 2; 3 4; 5 6];           % 3 time samples x 2 elements
d  = X(:);                      % stacked element by element
Dp = diag([1, -1i]);            % phase screen: element 2 delayed by a quarter period (phase -pi/2)
K  = kron(Dp, eye(3));          % 6x6, block diagonal
dab = K*d;
fprintf('kron(D, I3)*d = %s\n', mat2str(dab.', 3));
fprintf('equals vec(X*D): %d; correction kron(D'', I3) restores d: max err %.1e\n\n', ...
        isequal(dab, reshape(X*Dp, [], 1)), max(abs(kron(Dp', eye(3))*dab - d)));

%% 11. A tiny DAS operator: 2 pixels x (4 samples x 3 elements), linear interpolation
Nt = 4; Ne = 3;
tau = [2.5 2.0 2.5;             % pixel 1: delay per element, in samples (1-based, illustrative)
       3.25 3.0 3.25];          % pixel 2 (deeper)
rr = []; cc = []; ww = [];
for p = 1:2
    for e = 1:Ne
        n0 = floor(tau(p,e)); fr = tau(p,e) - n0;
        rr = [rr p p]; cc = [cc (e-1)*Nt+n0, (e-1)*Nt+n0+1]; ww = [ww 1-fr fr];
    end
end
keep = ww > 0;                  % drop zero weights (integer delays)
Bd = sparse(rr(keep), cc(keep), ww(keep), 2, Nt*Ne);
pm('full(B)', full(Bd));
fprintf('nnz(B) = %d of %d entries\n', nnz(Bd), numel(Bd));
dd = full(Bd(1,:)).';           % channel data of a unit scatterer at pixel 1 (H = B'' model)
pm('X = reshape(d, 4, 3)', reshape(dd, Nt, Ne));
fprintf('img = B*d = [%g %g]\n', full(Bd*dd));
pm('B*B'' (PSF matrix: column 1 = image of pixel 1)', full(Bd*Bd'));

%% 12. Why backslash and not inv: an ill-conditioned example
Hh = hilb(10); xt = ones(10,1); bh = Hh*xt;
xb = Hh \ bh; xi = inv(Hh)*bh;
fprintf('hilb(10): cond = %.1e\n', cond(Hh));
fprintf('  A\\b     : residual %.1e, error %.1e\n', norm(Hh*xb - bh)/norm(bh), norm(xb - xt)/norm(xt));
fprintf('  inv(A)*b: residual %.1e, error %.1e\n', norm(Hh*xi - bh)/norm(bh), norm(xi - xt)/norm(xt));

%% 13. Figure: structure of the operators
figure;
mats = {A, abs(R), real(F), T, abs(K), full(Bd)};
ttl  = {'A (general 3x3)', '|Z^H Z| (rank 2)', 'Re F_4', 'Toeplitz T (7x5)', '|kron(D, I_3)|', 'DAS B (2 x 12)'};
for k = 1:6
    subplot(2,3,k); imagesc(mats{k}); colorbar; title(ttl{k}); xlabel('column j'); ylabel('row i');
    if k < 6, axis image; end
end
colormap(gray);
