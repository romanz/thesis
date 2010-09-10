function beeper(f, t)
Fs = 8e3;
t = 0:1/Fs:t;
sound(sin(2*pi*f*t), Fs);