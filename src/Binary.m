% Abstract class for a binary operator
classdef Binary < Operator
properties
    op1, op2;
end

methods
    function self = Binary(op1, op2)
        if isa(op1, 'Operator') && isa(op2, 'Operator')
            assert(op1.grid == op2.grid);
        end
        if ~isa(op1, 'Operator')
            op1 = Const(op2.grid, op1);
        end
        if ~isa(op2, 'Operator')
            op2 = Const(op1.grid, op2);
        end
        self = self@Operator(op1.grid);
        self.op1 = op1;
        self.op2 = op2;
    end
end

end