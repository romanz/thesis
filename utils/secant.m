% SOL = SECANT(F, A, B, ITERS)
% Secant method to f(x)=0 starting at [a,b]. 
function sol = secant(f, x, iters)
    x1 = x(1); fprintf('%f -> ', x1); f1 = f(x1); fprintf('%e\n', f1);
    x2 = x(2); fprintf('%f -> ', x2); f2 = f(x2); fprintf('%e\n', f2);
    x3 = next(x1, x2, f1, f2);
    for i = 1:iters
        if x2 == x3 || abs(x2 - x3) < 1e-8
            break
        end
        fprintf('%f -> ', x3);
        f3 = f(x3);
        fprintf('%e\n', f3);
        x1 = x2; f1 = f2;
        x2 = x3; f2 = f3;
        x3 = next(x1, x2, f1, f2);
    end
    sol = x3;
end

function sol = next(x1, x2, f1, f2)
    dx = x2 - x1;
    dy = f2 - f1;
    if dx && dy
        sol = x1 - f1 * dx / dy;
    else
        sol = x2;
    end
end