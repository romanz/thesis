e = 0.1; 
n = 1e3; 
I = 0:n; 
h = 1/n;
d = h/(2*e); 
g = (1-d)/(1+d); 
u = (1 - g.^I) ./ (1 - g.^n); 
x = I/n; 
f = (exp(-x/e) - 1)/(exp(-1/e) - 1);
subplot 121; plot(x, u - f, 'Color', [0 0.6 0]); title('Error')
subplot 122; plot(x, u, 'b', x, f, 'r-.'); title('Solution')
