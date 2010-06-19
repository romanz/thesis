function I = checkerboard(sz)
    I = false(sz);
    I(1:2:end, 1:2:end) = true;
    I(2:2:end, 2:2:end) = true;
end
