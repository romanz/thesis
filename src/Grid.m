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
    x = x(:);
    xc = [ghost(flipud(x(1:3))); x; ghost(x(end-2:end))];
    xc = (xc(1:end-1) + xc(2:end)) / 2;
    switch type
        case 'central'
            x = xc;
        case 'interior'
            x = xc(2:end-1);
    end
end

function y = ghost(x)
    dx = diff(x);
    t = dx(2:end) ./ dx(1:end-1);    
    y = x(end) + dx(end) * t(end);    
%     y = x(end) + dx(end);
end
