function laplace(M ,T)
Nx = M;
Ny = M;
% sz = [Nx Ny];
% N = prod(sz);

ind = @(ix, iy) sub2ind([Nx Ny], ix, iy);

% S = stencil(:).';
% Dx = repmat([-1; 0; 1], [1 3]);
% Dy = Dx';
% Dx = Dx(:);
% Dy = Dy(:);
% 
% [Ix, Iy] = meshgrid(2:Nx-1, 2:Ny-1);
% Ix = Ix(:);
% Iy = Iy(:);
% Kx = repmat(Ix, [1 numel(Dx)]) + repmat(Dx', [numel(Ix) 1]);
% Ky = repmat(Iy, [1 numel(Dy)]) + repmat(Dy', [numel(Iy) 1]);
% K = ind(Kx, Ky);
% I = repmat(ind(Ix, Iy), [1 numel(stencil)]);
% Q = sparse(N, N);
% J = sub2ind([N N], I, K);
% S = repmat(S, [size(I, 1) 1]);
% Q(J) = S;
% L = sparse(N, N);
% for ix = 2:Nx-1
%     for iy = 2:Ny-1
%         J = ind(ix + Dx, iy + Dy);
%         L( ind(ix,iy),  J )  = stencil(:);
%     end
% end
% L-Q
% 
% R = sparse(N, N);
% Dx = repmat([-1; 0; 1], [1 3]);
% Dy = Dx';
% for i = 1:size(stencil, 1)
%     for j = 1:size(stencil, 2)
%         Kx = Ix + Dx(i, j);
%         Ky = Iy + Dy(i, j);
%         K = ind(Kx, Ky);
%         I = ind(Ix, Iy);
%         s = stencil(i, j);
%         R(sub2ind([N N], I, K)) = s;
%     end
% end
% L=R;
%% Try a simple way
s = [0 1.1 0; 1.2 -4 1.3; 0 1.4 0];
[Ix, Iy] = ndgrid(2:Nx-1, 2:Ny-1); % and not meshgrid!
Dx = repmat([-1; 0; 1], [1 3]);
Dy = Dx';
Dx = Dx(:);
Dy = Dy(:);
Ix = Ix(:);
Iy = Iy(:);

randn('state', 0);
U0 = round(randn(Nx, Ny));

[S, Kx, Ky] = stencil(s, Dx, Dy, Ix, Iy);
K = ind(Kx, Ky);
K = reshape(K, [[Nx Ny] - size(s) + 1, numel(s)]);
S = reshape(S, size(K));

V = U0;
tic;
for t = 1:T
    U = V(K);
    U = sum(S .* U, 3);
    V(2:end-1, 2:end-1) = U;
end
toc;
W = U0;
tic;
for t = 1:T
    W(2:end-1, 2:end-1) = conv2(W, fliplr(flipud(s)), 'valid');
end
toc;
norm(V - W, inf)
