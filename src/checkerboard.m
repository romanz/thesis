% I = CHECKERBOARD(SZ)
%   Create True/False checkerboard pattern of specified size.
function I = checkerboard(sz)
    I = false(sz);
    I(1:2:end, 1:2:end) = true;
    I(2:2:end, 2:2:end) = true;
end
