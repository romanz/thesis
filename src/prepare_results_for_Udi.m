function prepare_results_for_Udi

L = glob('MATs/513x513_Rmax=10/*beta=*_[513x513]_Rmax=10.0*0.mat');
S = cell(size(L));
i = 0;
for f = L
    i = i + 1;
    f = f{1};
    fprintf('.')
    S{i} = extract(load(f), 1);
    pause(0);
end

F = fieldnames(S{1});
N = numel(S);

res = struct();
for f = F'
    f = f{1};
    w = [];
    for i = 1:N
        s = S{i};
        v = s.(f);
        w = [w, v(:)];
    end
    res.(f) = w;
end

[~, I] = sort(res.beta);

for f = fieldnames(res)'
    f = f{1};
    w = res.(f);
    w = w(:, I);
    res.(f) = w;
end

for f = {'r', 't', 'zeta', 'alpha', 'Du'}
    f = f{1};
    w = res.(f);
    w = w(:, 1);
    res.(f) = w;
end
save('512x512_Rmax=10_zeta=10_alpha=0.5_Du=1_r=1.mat', 'res');

function res = extract(S, r)
res.r = r;
res.t = S.sol.grid.t;

g = @(op) regrid(Interp(Grid(res.r, res.t), op));
res.C = g(S.sol.C);
res.Phi = g(S.sol.Phi);
res.Vr = g(S.sol.Vr);
res.Vt = g(S.sol.Vt);

Er = -Deriv(Grid(S.sol.grid.r, S.sol.Phi.grid.t), S.sol.Phi, 1);
Et = -Deriv(Grid(S.sol.Phi.grid.r, S.sol.grid.t), S.sol.Phi, 2);

res.Er = g(Er);
res.Et = g(Et);

res.beta = S.sol.beta;
res.zeta = S.sol.zeta;
res.alpha = S.sol.alpha;
res.Du = S.sol.Du;
