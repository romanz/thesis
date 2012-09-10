function main_run

    N = 2.^(8);
    Rmax = 100;
    betas = 10.^[-2:0.25:0];
    A = 1;
    for a=A
        for n=N 
            for r=Rmax
                g = grids(r, n+1, n+1, a);
                main_ionex([], g, betas);
            end
        end
    end

end

% Create problem grids.
function [g] = grids(Rmax, Nr, Nt, a)
    r1 = logspace(0, log10(Rmax), Nr);
    r = Rgrid(Rmax, Nr, a);
    dr = diff(r(1:2));
    dt = dr;
    if nargin < 3
        Nt = ceil(pi / dt) + 1;
        Nt = Nt + ~mod(Nt, 2);
    end
    t = linspace(0, pi, Nt);
    r = r(:);
    t = t(:);
    rg = [2*r(1) - r(2); r; 2*r(end) - r(end-1)];
    tg = [2*t(1) - t(2); t; 2*t(end) - t(end-1)];
    rc = average(rg, [1; 1]/2);
    tc = average(tg, [1; 1]/2);
    
    g.Phi = Grid(rc, tc);
    g.C = Grid(rc, tc);
    g.Vr = Grid(r, tc);
    g.Vt = Grid(rc, t);
    g.P = Grid(rc(2:end-1), tc(2:end-1));
    g.r = r;
    g.t = t;
    
    dr = diff(r(1:2));
    dt = diff(t(1:2));
    fprintf('Nr:Nt = %d x %d\n', Nr, Nt)
    fprintf('dr:dt = %.4f:%.4f\n', dr, (mean(r(1:2) * dt)))
end

function r = Rgrid(Rmax, N, g)
    t = linspace(0, 1, N);
    if nargin >= 3
        t = t .^ g;
    end
    r = Rmax .^ t;
end
