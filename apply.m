function v = apply(v, R, g, T)
for t = 1:T
    v = R*v + g;
end
