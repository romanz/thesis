load
v = []; for k = 1:numel(solutions), v(k) = solutions{k}.Vinf; end
b = []; for k = 1:numel(solutions), b(k) = solutions{k}.beta; end
loglog(b, v)