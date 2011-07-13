function [res] = samegrid(g1, g2)
res = isequal(g1.r, g2.r) && isequal(g1.t, g2.t);
end