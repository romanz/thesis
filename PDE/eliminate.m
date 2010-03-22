function [A, f] = eliminate(A, f, J)
% Eliminate J variables from Av = f (for boundary elimination).
Af = [A f];
for k = 1:numel(J)
    j = J(k);
    I = find(Af(:, j));   
    % Apply Gaussian Elimination step (using an outer product)
    row = Af(j, :) / Af(j, j);
    col = Af(I, j);
    Af(I, :) = Af(I, :) - col * row;
end
A = Af(:, 1:end-1);
f = Af(:, end);
