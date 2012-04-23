function res = col(varargin)
n = 0;
for k = 1:numel(varargin)
    n = n + numel(varargin{k});
end
res = zeros(n, 1);
offset = 0;
for k = 1:numel(varargin)
    x = varargin{k};
    res(offset + (1:numel(x))) = x(:);
    offset = offset + numel(x);
end
