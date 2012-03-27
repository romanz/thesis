function value = regrid(op)
    value = reshape(op.res(), op.grid.size);
end
