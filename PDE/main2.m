function main2
    randn('state', 1);
    rand('twister', 1);
    fprintf([repmat('-', 1, 80) '\n']); 
    sz = 1+2.^([1 1]*3);
    N = prod(sz);
    
    %% Create the grid
    x = linspace(-1, 1, sz(1));
    y = linspace(-1, 1, sz(2));
    [X, Y] = ndgrid(x, y);
    I = interior(sz);
    Bl = boundary(I, [-1 0]); Il = shift(Bl, [+1 0]);
    Br = boundary(I, [+1 0]); Ir = shift(Br, [-1 0]);
    Bd = boundary(I, [0 -1]); Id = shift(Bd, [0 +1]);
    Bu = boundary(I, [0 +1]); Iu = shift(Bu, [0 -1]);

    Vx = stagger(I, X, Y, 1, @(x, y) y);
    Vy = stagger(I, X, Y, 2, @(x, y) x);

    %% Solve the coupled problem
    U = randn(sz);
    C = rand(sz);
    
    Jcorners = sub2ind(sz, [1 1 sz(1) sz(1)], [1 sz(2) 1 sz(2)]);
    for t = 1:300
        [A] = laplacian(I, X, Y);
        f = zeros(sz); % Right-Hand Side
        [A, f] = dirichlet(A, f, find(Bl), exp(-U(Il)));
        [A, f] = dirichlet(A, f, find(Br), 1 + 0*Y(Br));
        [A, f] = neumann(A, f, find(Bd), find(Id), 0*X(Bd));
        [A, f] = neumann(A, f, find(Bu), find(Iu), 0*X(Bu));
        C = solve(A, f, I, C);

        C(Jcorners) = NaN; subplot 121; mesh(X, Y, C); title('C')

        [A] = laplacian(I, X, Y, C) - 1*advection(I, X, Y, Vx, Vy, 'upwind');
        f = zeros(sz); % Right-Hand Side
        [A, f] = dirichlet(A, f, find(Bl), -log(C(Il)));
        [A, f] = dirichlet(A, f, find(Br), 0*Y(Br));
        [A, f] = neumann(A, f, find(Bd), find(Id), 0*X(Bd));
        [A, f] = neumann(A, f, find(Bu), find(Iu), 0*X(Bu));
        U = solve(A, f, I, U);        

        U(Jcorners) = NaN; subplot 122; mesh(X, Y, U); title('\Phi')
    end
%     reshape(C(I), sz-2), reshape(U(I), sz-2)
    reshape(C, sz) - 1, reshape(U, sz)

end

function v = solve(A, f, I, v0)
    N = numel(f);
    f = f(:);
    [Ai, fi] = eliminate(A, f, find(~I));
    [Ai, fi] = restrict(Ai, fi, I);
    v = zeros(size(f));
    [R, T, d] = jacobi(Ai, fi);
    v(I) = iterate(v0(I), T, d, 100, 'Jacobi', 'silent');
%     v(I) = Ai \ fi;
    
    A(I, :) = 0;
    A(sub2ind([N N], find(I), find(I))) = 1;
    f(I) = v(I);
    
    [A, f] = eliminate(A, f, find(I));
    assert(nnz(diag(diag(A)) - A) == 0);
    J = find(diag(A)); 
    % Verify that boundary coefficient is 1.
    assert(all(diag(A(J, J) == 1)))
    % Update boundary values for the solution.
    v(J) = f(J); 
    v = reshape(v, size(I));
end
