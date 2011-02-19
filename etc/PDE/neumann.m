function [A, f] = neumann(A, f, J, I, dU)
for k = 1:numel(J)
    j = J(k);
    A(j, [j I(k)]) = [1, -1];
    f(j) = dU(k);
end
