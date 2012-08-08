clear;

U1 = 4.25752957407993;
U3 = 0.966869;

S{1} = load('20120805125313');
S{2} = load('20120805125331');
S{3} = load('20120805125512');

Vn = zeros(numel(S), 1);
for k = 1:numel(S)
    Vn(k) = S{k}.V;
    g = S{k}.g;
end
Vr = Vn(3)*4/3 - Vn(2)/3;

b = S{3}.betas;
V1 = U1*b;
V3 = U3*b.^3;

plot([b b b], Vn, 'o', b, Vr, 's', b, V1, 'd', b, V1+V3, 'd')
legend('Numerical', 'Richardson', 'Linear', 'Cubic')
dVn = diff(Vn)
dVn(1) / dVn(2)
