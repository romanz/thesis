function M = redblack(sz, A)
[I, J] = ndgrid(1:sz(1), 1:sz(2));
K = mod(I + J, 2);
I0 = num2cell(find(K == 0));
I1 = num2cell(find(K == 1));
M{1} = vanka(A, I0, I0);
M{2} = vanka(A, I1, I1);
