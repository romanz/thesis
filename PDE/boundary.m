function B = boundary(I, P)

B = ~I & circshift(I, P);
