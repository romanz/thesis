function main2
    randn('state', 0);
    fprintf([repmat('-', 1, 80) '\n']); 
    sz = 1+2.^([1 1]*2);
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

    %% Solve the coupled problem
    U = 0 + 0.1*(randn(sz));
    for t = 1:100
        [A] = laplacian(I, X, Y);
        f = zeros(sz); % Right-Hand Side
        Cl = exp(-U(2:end-1, 1));
        Cr = ones(sz(1)-2, 1);
        [A, f] = dirichlet(A, f, find(Bl), Cl);
        [A, f] = dirichlet(A, f, find(Br), Cr);
        [A, f] = neumann(A, f, find(Bu), find(Iu), zeros(nnz(Bu)));
        [A, f] = neumann(A, f, find(Bd), find(Id), zeros(nnz(Bd)));
        C = solve(A, f, I);

        [A] = laplacian(I, X, Y, C);
        f = zeros(sz); % Right-Hand Side
        Ul = -log(C(2:end-1, 1));
        [A, f] = dirichlet(A, f, find(Bl), Ul);
        [A, f] = neumann(A, f, find(Br), find(Ir), zeros(nnz(Br)));
        [A, f] = neumann(A, f, find(Bu), find(Iu), zeros(nnz(Bu)));
        [A, f] = neumann(A, f, find(Bd), find(Id), zeros(nnz(Bd)));
        U = solve(A, f, I);        
    end
    reshape(C, sz), reshape(U, sz)

end

function v = solve(A, f, I)
    N = numel(f);
    f = f(:);
    [Ai, fi] = eliminate(A, f, find(~I));
    [Ai, fi] = restrict(Ai, fi, I);
    v = zeros(size(f));
    v(I) = Ai \ fi;
    
    A(I, :) = 0;
    A(sub2ind([N N], find(I), find(I))) = 1;
    f(I) = v(I);
    
    [A, f] = eliminate(A, f, find(I));
    assert(nnz(diag(diag(A)) - A) == 0);
    J = find(diag(A));
    assert(all(diag(A(J, J) == 1)))
    v(J) = f(J);
end
