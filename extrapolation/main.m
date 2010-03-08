% Main simulation script
%
% Roman Zeyde, Computer Science Department
% Technion -- Israel Institute of Technology
% romanz@cs.technion.ac.il

N = 1000; % Problem dimension
T = sparse(N, N);
T(sub2ind([N N], 1:N, 1:N))   = [5 repmat(6, [1 N-2]) 5];
T(sub2ind([N N], 2:N, 1:N-1)) = [2 repmat(3, [1 N-3]) 2];
T(sub2ind([N N], 3:N, 1:N-2)) = 1;
T(sub2ind([N N], 4:N, 1:N-3)) = 1;
T(sub2ind([N N], 1:N-1, 2:N)) = [2 repmat(3, [1 N-3]) 2];
T(sub2ind([N N], 1:N-2, 3:N)) = 1;
T(sub2ind([N N], 1:N-3, 4:N)) = 1;
T = 0.06 * T; % T is NxN septa-diagonal SPD matrix.
sol = ones(N, 1); % Pre-defined solution
d = (eye(N) - T) * sol; % sol = T * sol + d
k = 10; % Extrapolation order
L = 8; % # of cycles
w = 2; % Smoothing factor
% Iterative method: 
A = w*T + (1-w)*speye(N);
b = w*d;
F = @(x) A*x + b; % w(T*x + d) + (1-w)x
s = zeros(N, 1); % Initial guess
for i = 1:20
    s = F(s);
end
err = @(x) norm(sol - x);
[t, r1] = extrapolate(s, F, k, L, 'MPE');
[t, r2] = extrapolate(s, F, k, L, 'RRE');
