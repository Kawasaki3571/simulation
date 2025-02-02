p0 = 1;
r = 1;
cal_time = 0.18;
dt = cal_time / 90000;
dx = c/(2000*100);
t = 0 : dt : cal_time - dt;
freq_start = 1000;
freq_add = 3000;
c = 340;
pai = 3.1415;
sweep = 1;
gosahyoka = 0;

if sweep == 1
    freq = freq_start + freq_add*t/0.18;
else
    freq_start = 2000;
    freq_add = 0;
    freq = freq_start + freq_add*t/0.18;
end

omega = freq * 2 * pai;
hasu = 2 * pai * freq / c;
j = sqrt(-1);
xL = 0.5;
x = 0 : 0.001 : 2;

for t_x = 1 : cal_time/dt
     time = t_x * dt;
    x = 0 : dx : (time*c);
    p = p0.*exp(j .* (hasu(t_x) .*x -omega(t_x) .* time ) );
    if mod(t_x, 50) == 0
        plot(x, p)
        drawnow
        pause(0.15)
    end
%    plot(x, p)
%    drawnow
end
