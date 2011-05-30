function r = Rgrid(d, q, Rmax)
r = 1;
while r(end) < Rmax
    r(end+1) = r(end) + d;
    d = d * q;
end
