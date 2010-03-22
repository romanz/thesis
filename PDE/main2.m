function main2
    randn('state', 0);
    fprintf([repmat('-', 1, 80) '\n']); 
    sz = 1+2.^([1 1]*2);
    N = prod(sz);
    %% Create the grid
    x = linspace(-1, 1, sz(1));
    y = linspace(-1, 1, sz(2));
    [X, Y] = ndgrid(x, y);

    %% Solve the coupled problem
    U = 0 + 0.1*(randn(sz));
    for t = 1:100
        [A, f, I] = C_prob(sz, X, Y, exp(-U(2:end-1, 1)), ones(sz(1)-2, 1));
        C = solve(A, f, I);
        [A, f, I] = U_prob(sz, X, Y, C, -log(C(2:end-1, 1)));
        U = solve(A, f, I);        
    end
    reshape(C, sz), reshape(U, sz)

    function [A, f, I] = C_prob(sz, X, Y, Cl, Cr)
        [L, M, I] = laplacian(sz, X, Y);
        Bl = boundary(I, [-1 0]); Il = shift(Bl, [+1 0]);
        Br = boundary(I, [+1 0]); Ir = shift(Br, [-1 0]);
        Bd = boundary(I, [0 -1]); Id = shift(Bd, [0 +1]);
        Bu = boundary(I, [0 +1]); Iu = shift(Bu, [0 -1]);
        A = dinv(M) * L;
        f = zeros(sz); % Right-Hand Side

        [A, f] = dirichlet(A, f, find(Bl), Cl);
        [A, f] = dirichlet(A, f, find(Br), Cr);
        [A, f] = neumann(A, f, find(Bu), find(Iu), zeros(nnz(Bu)));
        [A, f] = neumann(A, f, find(Bd), find(Id), zeros(nnz(Bd)));
    end

    function [A, f, I] = U_prob(sz, X, Y, C, Ul)
        [L, M, I] = laplacian(sz, X, Y, C);
        Bl = boundary(I, [-1 0]); Il = shift(Bl, [+1 0]);
        Br = boundary(I, [+1 0]); Ir = shift(Br, [-1 0]);
        Bd = boundary(I, [0 -1]); Id = shift(Bd, [0 +1]);
        Bu = boundary(I, [0 +1]); Iu = shift(Bu, [0 -1]);
        A = dinv(M) * L;
        f = zeros(sz); % Right-Hand Side

        [A, f] = dirichlet(A, f, find(Bl), Ul);
        [A, f] = neumann(A, f, find(Br), find(Ir), zeros(nnz(Br)));
        [A, f] = neumann(A, f, find(Bu), find(Iu), zeros(nnz(Bu)));
        [A, f] = neumann(A, f, find(Bd), find(Id), zeros(nnz(Bd)));
    end

    function v = solve(A, f, I)
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
end
