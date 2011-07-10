function C = dot_prod(A, B)
C = A * B;
if isequal(A, 0) || isequal(B, 0)
    C = sparse(0);  
end
