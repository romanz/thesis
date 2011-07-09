classdef Operator < handle
properties
    grid; % output
end
methods
    function op = Operator(grid)
        op.grid = grid;
    end
    function op = plus(op1, op2)
        op = Sum(op1, op2);
    end
    function op = minus(op1, op2)
        op = Sum(op1, -op2);
    end
    function op = uminus(op)
        op = Product(-1, op);
    end
    function op = mtimes(op1, op2)
        op = Product(op1, op2);
    end
    function op = log(op)
        op = Func(op.grid, op, 'log(x)');
    end
    function op = exp(op)
        op = Func(op.grid, op, 'exp(x)');
    end
    function op = sinh(op)
        op = Func(op.grid, op, 'sinh(x)');
    end
    function op = cosh(op)
        op = Func(op.grid, op, 'cosh(x)');
    end
end
end