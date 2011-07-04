function g = Grid(radius, theta)
    g.r = radius(:);
    g.t = theta(:);
    [g.R, g.T] = ndgrid(g.r, g.t);
    g.sz = [numel(g.r), numel(g.t)];
    g.numel = prod(g.sz);
    function res = init(varargin)
        res = Grid(axis(g.r, varargin{1}), axis(g.t, varargin{2}));
        res.ghost = ~strcmp('interior', varargin);
    end
    g.init = @init;
end

function x = axis(x, type)    
    if isempty(type)
        return
    end
    xc = [[2,-1]*x(1:2); x; [-1,2]*x(end-1:end)];
    xc = (xc(1:end-1) + xc(2:end)) / 2;
    switch type
        case 'central'
            x = xc;
        case 'interior'
            x = xc(2:end-1);
    end
    disp(x)
end
