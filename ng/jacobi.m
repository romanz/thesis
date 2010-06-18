function M = jacobi(sz, A)
I = num2cell(find(true(sz)));
M{1} = vanka(A, I, I);
