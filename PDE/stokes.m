% 2D Stokes equation discretization.
function [A, f] = stokes(sz, X, Y, Fx, Fy, VxL, VxR, VyD, VyU)

I = interior(sz);
Gx_P = gradient(sz - 2, X(I), Y(I), 1);
Gy_P = gradient(sz - 2, X(I), Y(I), 2);
Ni = nnz(I); % # of interior points
    
% X-staggered grid
h = [1; 1]/2;
Xs = average(X, h);
Ys = average(Y, h);
I = interior(sz - [1 0]);
L_Vx = laplacian(I, Xs, Ys);
L_Vx = L_Vx(I, :);
Gx_Vx = gradient(sz - [1 0], Xs, Ys, 1);
I = interior(sz - [2 0], [0 1]);
Gx_Vx = Gx_Vx(I, :);

% Y-staggered grid
h = [1, 1]/2;
Xs = average(X, h);
Ys = average(Y, h);
I = interior(sz - [0 1]);
L_Vy = laplacian(I, Xs, Ys);
L_Vy = L_Vy(I, :);
Gy_Vy = gradient(sz - [0 1], Xs, Ys, 2);
I = interior(sz - [0 2], [1 0]);
Gy_Vy = Gy_Vy(I, :);

% - lapl V[x] + grad[x] P = F[x]
% - lapl V[y] + grad[y] P = F[y]
% grad[x] V[x] + grad[y] V[y] = 0
% u = [Vx; Vy; P];
A = [[-L_Vx, sparse(size(L_Vx, 1), size(L_Vy, 2)), Gx_P]; ...
     [sparse(size(L_Vy, 1), size(L_Vx, 2)), -L_Vy, Gy_P]; ...
     [Gx_Vx, Gy_Vy, sparse(Ni, Ni)]];
