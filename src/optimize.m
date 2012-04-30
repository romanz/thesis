% Operator expression optimizer
function op = optimize(op)
%     disp(class(op))
    if isa(op, 'Join')
        for k = 1:numel(op.ops)
            op.ops{k} = optimize(op.ops{k});
        end
        return
    end
    if isa(op, 'Binary')
        op.op1 = optimize(op.op1);
        op.op2 = optimize(op.op2);
        if isa(op, 'Product')
            if isa(op.op1, 'Linear') && isa(op.op2, 'Const')
                op = Linear(op.grid, op.op1.op, spdiag(op.op2.val) * op.op1.L);
                return
            end
            if isa(op.op1, 'Const') && isa(op.op2, 'Linear')
                op = Linear(op.grid, op.op2.op, spdiag(op.op1.val) * op.op2.L);
                return
            end
        end
        return
    end
    if isa(op, 'Linear')
        op.op = optimize(op.op);
        if isa(op.op, 'Linear')
            op.grid = op.op.grid;
            op.L = op.L * op.op.L;
            op.op = op.op.op;
        end
        return
    end
    if isa(op, 'Func')
        op.op = optimize(op.op);
        return
    end    
end

