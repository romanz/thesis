function [J, I] = boundary(interior, P)
J = ~interior & shift(interior, P);
I = shift(J, -P);
