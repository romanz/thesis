function value = regrid(op)
value = reshape(op.res(), op.grid.sz);
end