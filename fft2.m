r = 0.5;
p0 = 1;
xL = 0.5;
c = 340;
Fs = 2000;
Ts = 1 / Fs;
cal_time = 0.18;
L = cal_time / Ts;
k = 1;
j = sqrt(-1);
x = 0 : 0.0001 : 0.1;
xn = 2 / 0.001 + 1;
dt = 0.0001;
t = 0 : dt : cal_time;


i = 1 : 17888;
M = csvread('二万点データ.csv',2,1);
diff = M(i,4);
% disp(diff)

% Q = dt * cumtrapz(pin);
freq = 2000 + 7000 * i /17888;

disp(freq(1))
disp("あああ")
disp(size(freq))
% disp(size(spec))
disp(size(x))
df = 0.5;
w = 2 * pi * freq;
spec = p0^2 *(1 + r^2 + 2 * r * cos(2*pi*(2*xL*freq/c)));

disp("A")
spec_s = {};


% plot(freq,exp(-1*j*2*pi*freq))
for f = 1 : 17888
%     for n = 1 : xn
%         spec_s(n,f) = spec(f) .* exp(-1*j*x(n).*freq(f));
%         fun = @(x) spec(f) .* exp(-1*j*x.*freq(f));
%         spec_s{f} = fun;
%     end
%     diff(f) = diff(f).*exp(-1*j*x.*freq(f));
    fun = @(x) diff(f) .* exp(-1*j*x.*freq(f));
    spec_s{f} = fun;
end
disp(size(freq))

% plot(x,spec_s{1000}(x))

disp("B")
disp(size(diff))


fen = @(x) 0 ;
for f = 1 : 17888
%     fen = @(x) fen(x) + df * spec_s{f}(x);
    fen = @(x) fen(x) + df * spec_s{f}(x);
%       fen = @(x) df * cumtrapz(spec_s{f}(x));
    f
    if f == 17800
%         x = 0 : 0.001 : 2;
%         plot(x, fen(x*4*pi/c))
        break;
    end
end
% fen(10)
        x = 0 : 0.001 : 2;
        plot(x, fen(x*4*pi/c))

% plot(x, fin)
x = 0 :0.001: 1;
% plot(x,fin(x))
% fin(10)
% plot(x, fin(x))

disp("owari")


Y = fft(spec);
P2 = abs(Y/L);
P1 = P2(1:L/2 + 1);
p = 1 : L/2 + 1;
f = Fs*(0:(L/2))/L;
% plot(x, P)