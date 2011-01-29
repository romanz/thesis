% Splits x into seperate variables according to specified sizes.
%   [x1, x2, x3] = split(x, sz1, sz2, sz3);
%   where sz{i} = size(x{i}) and x = [x1(:); x2(:); x3(:)];
function [varargout] = split(x, varargin)
N = numel(varargin);
varargout = cell(N);
offset = 0;
for k = 1:N
    sz = varargin{k};
    m = prod(sz);
    varargout{k} = reshape(x(offset + (1:m)), sz);
    offset = offset + m;
end
assert(offset == numel(x));
