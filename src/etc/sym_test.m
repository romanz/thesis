clear
syms r t
Dt = @(f) diff(f, t);
Dr = @(f) diff(f, r);
Vr = cos(t)/r^3;
Vt = 0.5*sin(t)/r^3;
r^(-2)*Dr(r^2*Dr(Vt))
r^(-2)*Dt((1/sin(t))*Dt(Vt*sin(t)) + 2*Vr)