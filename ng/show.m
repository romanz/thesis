function show(M, msg)
% M(1, 1) = NaN;
% M(1, end) = NaN;
% M(end, 1) = NaN;
% M(end, end) = NaN;
mesh(M.')
title(msg)
xlabel('R')
ylabel('\theta')
