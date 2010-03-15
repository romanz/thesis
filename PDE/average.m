function A = average(A, h)
A = convn(A, h, 'valid');
