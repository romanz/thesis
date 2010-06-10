function A = shave(A, rows, cols)
A = A(1+rows : end-rows, 1+cols : end-cols);
