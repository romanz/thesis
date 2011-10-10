function L = vlapl(F)
    syms r t;
    Dr = @(f) diff(f, r);
    Dt = @(f) diff(f, t);
    sint = sin(t);
    Vr = F(1);
    Vt = F(2);
    L = [Dr(r^(-2) * Dr(r^2 * Vr)) + (1/(r^2 * sint)) * Dt((Dt(Vr) - 2*Vt) * sint);
         r^(-2) * Dr(r^2 * Dr(Vt)) + r^(-2) * Dt((1/(sint))*Dt(Vt * sint) + 2*Vr)];
    L = simple(L);
end
