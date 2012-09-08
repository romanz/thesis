function S = load_results(s)

p = sprintf('%s/', s);
D = dir([p '*x*.mat']);

if isempty(D)
    fprintf('No MAT files found.\n')
    return
end

S = struct;
for k = 1:numel(D)
    d = D(k);
    f = d.name;
    s = load([p, f], 'betas', 'k', 'v', 'g');
    disp([p f])
    S(k).beta = s.betas(s.k);
    S(k).v = s.v(end);
    S(k).n = numel(s.g.r)-1;
end

b = unique([S.beta]);
n = unique([S.n]);
R = {};
for k = 1:numel(S)
    i = find(b == S(k).beta);
    j = find(n == S(k).n);
    R{i,j} = S(k);
end

Vn = cellfun(@(r) r.v, R);
b = b(:);
n = n(:);
Vr = Vn(:, end-1:end) * [-1; 4]/3;
save([p 'data.mat'], 'Vn', 'b', 'n', 'Vr')