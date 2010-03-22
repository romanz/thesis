function B = boundary(I, P)
B = ~I & shift(I, P);
