% Unit Testing for spinit

%% Simple case
n = 5; 
I = [1:n; 1:n; 1:n]; 
J = [1:n; 2:n+1; 3:n+2];
mask = spinit(I, J, [n n+2]);
for e = [0 1 2]
    values = repmat([1+e;-2; 1-e], [1 n]);
    S = mask(values);
    q = (S == sparse(I, J, values, n, n+2));
    assert(all(q(:)))
end

%% Large matrix
n = 1e3;
[I, J] = find(sprand(n, n, 5/n));
mask = spinit(I, J, [n n]);
for i = 1:10
    values = randn(size(I));
    S = mask(values);
    q = (S == sparse(I, J, values, n, n));
    assert(all(q(:)))
end

%% Compare speed up
n = 10e3;
[I, J] = find(sprand(n, n, 5/n));
L = numel(I);
x = 1:L;
T = 1e3;
tic;
mask = spinit(I, J, [n n]);
for i = 1:T
    values = x + i;
    S = mask(values);
end
t0 = toc;
tic;
for i = 1:T
    values = x + i;
    S = sparse(I, J, values, n, n);
end
t1 = toc;
fprintf('Speed-Up: %.1f%%\n', 100 * t1 / t0);