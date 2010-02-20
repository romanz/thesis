%% Test which stencil is faster
randn('state', 0);
N = 2^6+1;
V0 = sign(randn(N, N));
H1 = 0.0*randn(N-2, N-2, 1)+1/2;
H2 = 0.0*randn(N-2, N-2, 1)+1/4;
H3 = 0.0*randn(N-2, N-2, 1)+1/8;
H4 = 0.0*randn(N-2, N-2, 1)+1/16;
H5 = 0.0*randn(N-2, N-2, 1)+1/32;
T = 10;

% Block
V = V0;
tic
I = 2:N-1;
Ip= I+1;
Im= I-1;
for t = 1:T
    U = V(I,  I)  .* H1 + ...
        V(I,  Im) .* H2 + ...
        V(I,  Ip) .* H3 + ...
        V(Im, I)  .* H4 + ...
        V(Ip, I)  .* H5;
    V(I, I) = (U);
end
toc
V1 = V;

% Sparse
sz = [N N];
S = 5;
V = V0;
[I, J] = ndgrid(1:N-2, 1:N-2);
I = I(:); J = J(:);
P = sub2ind(sz, ...
    repmat(I, [1 S]) + repmat([-1 0 0 0 1], size(I)) + 1, ...
    repmat(J, [1 S]) + repmat([0 -1 0 1 0], size(J)) + 1);
K = sub2ind(sz, repmat(I, [1 S]) + 1, repmat(J, [1 S]) + 1);
H = [H4(:) H2(:) H1(:) H3(:) H5(:)];
A = sparse(K(:), P(:), H(:), N*N, N*N);
I = any(A, 2);
A = A(I, :);
V = V(:);
tic
for t = 1:T
    U = A * V;
    V(I) = (U);
end
toc
V2 = reshape(V, sz);

% Show diff.
E = V1 - V2;
imagesc(E); 
colorbar
norm(E(:)) / norm(V1(:))
