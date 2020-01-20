pi = 3.1416;
t = 0 : 0.01 : 10;
freq = 10000;
sina = sin(2*pi*freq*t);
sinb = sin(2*pi*freq*t + pi/3);
sinr = sina + sinb;
plot(t, sina)