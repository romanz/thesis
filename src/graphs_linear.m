clear;
clf;
hold on;
res(1) = show('Linear/pos/2012Mar19_192151.mat', '-');
res(end+1) = show('Linear/pos/2012Mar19_180958.mat', '--');
res(end+1) = show('Linear/pos/2012Mar19_173821.mat', '-.');
res(end+1) = show('Linear/pos/2012Mar19_175306.mat', ':');

legend([res.h], {res.s}, 'FontSize', 16, 'Location', 'SouthEast')
print -deps linearPos.eps


clear;
clf;
hold on;
res(1) = show('Linear/neg/2012Mar19_182449.mat', ':');
res(end+1) = show('Linear/neg/2012Mar19_183925.mat', '-.');
res(end+1) = show('Linear/neg/2012Mar19_172418.mat', '--');
res(end+1) = show('Linear/neg/2012Mar19_190517.mat', '-');

legend([res.h], {res.s}, 'FontSize', 16, 'Location', 'NorthEast')
print -deps linearNeg.eps

