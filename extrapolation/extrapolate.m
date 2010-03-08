% Extrapolate iterative method F using specified method.
% Use L cycles of k-order extrapolation, starting from x0.
%
% Roman Zeyde, Computer Science Department
% Technion -- Israel Institute of Technology
% romanz@cs.technion.ac.il

function [x0, residuals] = extrapolate(x0, F, k, L, method)
    N = numel(x0);
    Q = zeros(N, k+1);
    switch method
        case 'MPE',  method = @mpe;
        case 'RRE',  method = @rre;
        otherwise ,  error('Unknown method (%s)!', method);
    end
    % Perform L cycles of extrapolation method
    residuals = zeros(L, 1);
    for t = 1:L
        % Compute (k+2) vectors, including x0
        Q(:, 1) = F( x0 );        
        residuals(t) = norm(x0 - Q(:, 1), 2); 
        for i = 1:k 
            Q(:, i+1) = F( Q(:, i) );
        end
        % Compute differences (k+1)
        for i = k:-1:1 
            Q(:, i+1) = Q(:, i+1) - Q(:, i);
        end
        Q(:, 1) = Q(:, 1) - x0;
        % Perform QR decomposition
        [Q, R] = MGS(Q); 
        % Perform extrapolation
        [gamma] = method(R, k); % s.t. x0 = X * gamma
        xi = 1 - cumsum(gamma(1:k)); % s.t. x0' = x0 + U * xi
        eta = R(1:k, 1:k) * xi; % since U = Q * R
        x0 = x0 + Q(:, 1:k) * eta; % s.t. x0' = x0 + Q * R * xi
    end
end

% Minimal Polynomial Extrapolation
function [gamma, residual] = mpe(R, k)
    c = backsub(R(1:k, 1:k), -R(1:k, k+1));
    c = [c; 1];
    gamma = c / sum(c);
    residual = abs(gamma(end)) * R(end, end);
end

% Reduced Rank Extrapolation
function [gamma, residual] = rre(R, k)
    e = ones(k+1, 1);
    d = backsub(R, backsub(R', e));  
    lambda = 1 / sum(d);
    gamma = lambda * d;
    residual = sqrt(lambda);
end

