function tone(F, T)
Fs = 48e3;
t = 0:(1/Fs):T;
x = sin(2*pi*F*t) .* sin(pi*t/t(end)) * 0.5;
wavplay(x, Fs)
