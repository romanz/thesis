function I = shift(I, Q)

sz = size(I);
Qp = max(Q, 0); % positive shifts
Qn = min(Q, 0); % negative shifts

I = [zeroes(Q(1), size(I, 2)); ...
     I(1-Qn(1):end-Qp(1), :); ...
     zeroes(-Q(1), size(I, 2))];

I = [zeroes(size(I, 1), Q(2)), ...
     I(:, 1-Qn(2):end-Qp(2)), ...
     zeroes(size(I, 1), -Q(2))];

function Z = zeroes(n, m)
    Z = zeros(min(max([n m], 0), sz));
end

end