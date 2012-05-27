function [P] = upwind(gridV, V, dim)
    dir = (1:2 == dim);
    assert(any(dir));
    h = ones(dir + 1);
    q = ~dir; % remove other dimension's velocity ghost points
    I = gridV.I(1+q(1):end-q(1), 1+q(2):end-q(2));
    K1 = find(I | shift(I, -dir)); % positive flow face
    K2 = find(I | shift(I, +dir)); % negative flow face
    K = logical( convn(I, h, 'valid') ); % interior grid
    % The following code maybe optimized using SPINIT scheme.
    Y = convn(V, h, 'valid') >= 0;
    % Upwind logic:
    % > Positive flow, Y = 1, take K1 face.
    % < Negative flow, Y = 0, take K2 face.
    Y = Y(K);
    P = sparse([1:nnz(K), 1:nnz(K)], [K1, K2], [Y, ~Y], nnz(K), numel(I));
end
