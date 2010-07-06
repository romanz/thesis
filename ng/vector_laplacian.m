clear
fprintf('-------------------------------------------------------\n')

syms r t a
Vr = cos(t)*(1 - 3*a/(2*r) + a^3/(2*r^3));
Vt = -sin(t)*(1 - 3*a/(4*r) - a^3/(4*r^3));

subs(Vr, r, a)
subs(Vr, r, inf)

subs(Vt, r, a)
subs(Vt, r, inf)

Lr = diff(r^2 * diff(Vr, r), r) / r^2 + ...
     diff(sin(t) * diff(Vr, t), t) / (r^2 * sin(t)) - ...
     2 * (Vr + diff(sin(t) * Vt, t) / sin(t)) / r^2;
Lr = simple(Lr)

Lt = diff(r^2 * diff(Vt, r), r) / r^2 + ...
     diff(sin(t) * diff(Vt, t), t) / (r^2 * sin(t)) ...
     - Vt / (r * sin(t))^2 + 2*diff(Vr, t) / r^2;
Lt = simple(Lt)

P = -3*a/(2*r^2) * cos(t);

Fr = diff(P, r);
Fr = simple(Fr)

Ft = diff(P, t) / r;
Ft = simple(Ft)

(Lr == Fr && Lt == Ft)

div = diff(r^2 * Vr, r) / r^2 + diff(sin(t) * Vt, t) / (r * sin(t));
div = simple(div)

Srr = 2*diff(Vr, r) - P; Srr = subs(Srr, r, a)
Srt = diff(Vr, t) / r + r * diff(Vt/r, r); Srt = subs(Srt, r, a)
Fz = Srr*cos(t) - Srt*sin(t);
Fz = simple(Fz * 4*pi*a^2)