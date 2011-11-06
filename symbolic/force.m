function [f, fn, fm] = force(v, p, phi)
    syms r t pi
    Vr = v(1);    Vt = v(2);
    Fn = simple([-p + 2*Dr(Vr); Dr(Vt) + (Dt(Vr) - Vt)/r]);
    Fm = simple([((Dr(phi)^2 - (Dt(phi)/r)^2))/2; Dr(phi)*Dt(phi)/r]);
    
    axial = @(F) simple(subs([cos(t), -sin(t)] * F, r, 1));
    total = @(F) simple(int(axial(F) * 2*pi*sin(t), t, 0, pi));
    f = total(Fn + Fm);
    fn = total(Fn);
    fm = total(Fm);
end

function DfDr = Dr(f)
    DfDr = diff(f, 'r');
end

function DfDt = Dt(f)    
    DfDt = diff(f, 't');
end
