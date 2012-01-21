% C = dot_prod(A, B)
% If A or B are zero scalars, C is sparse zero scalar.
function C = dot_prod(A, B)

    if isequal(A, 0) || isequal(B, 0)
        C = sparse(0);  
    else
        C = A * B;
    end

end