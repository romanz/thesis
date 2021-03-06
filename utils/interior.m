function I = interior(sz, dims)

if nargin < 2
    dims = true(numel(sz), 1);
end

% Array of true values
I = true(sz);

% Exclude boundary points' indices (special handling of 1D)
if sz(1) > 1 && dims(1)
    I([1 sz(1)], :) = false;
end
if sz(2) > 1 && dims(2)
    I(:, [1 sz(2)]) = false; 
end
