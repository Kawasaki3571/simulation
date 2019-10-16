r = 0.5;
p0 = 1;
xL = 2;
c = 340;
L = 1500;
Fs = 1000;
Ts = 1 / Fs;
t = Ts*(1:L);
f = 1000 + t * 6 / Ts;
spec = p0^2*(1 + r^2 + 2*r*cos(2*pi*(2*xL*f)/c));
% plot(f, spec)
Y = fft(spec,L);
P2 = Y;
P1 = P2(1:L/2 + 1);
% disp(size(P1))
p = 1 : 751;
plot(p, P1)
disp(8 * pi / 340);