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