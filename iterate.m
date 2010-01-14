% Iterate T times using B on residual(v).
% [V, HISTORY] = ITERATE(V, B, RESIDUAL, T)
function [v, history] = iterate(v, B, residual, T)
history = zeros(numel(v), T);
for t = 1:T
    history(:, t) = v;
    r = residual(v);
    v = v + B * r;
end
