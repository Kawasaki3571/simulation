mode_plot = 3;
hekomi_bool = 1;
sweep = 1;

c = 340;
rou = 1.2;
l = 52;
x_in = 50;
hekomi = 51;
freq = 2000;
freq_start = 1000;
freq_add = 3000;
ramuda = c / freq ;
dx = ramuda / 100;
crn = 0.4;

dt = dx * crn / c;
ix = round(l / dx);
id = round(x_in / dx);
sokutei_end_x  = 52;
end_i = sokutei_end_x / dx;
end_point = round(end_i);
hekomi_point = round(hekomi/dx);
p1 = zeros(1,ix + 1);
p2 = zeros(1,ix + 1);
u1 = zeros(1, ix + 1);
u2 = zeros(1,ix + 1);
cal_time = 0.18;

time_max_g = round(cal_time/dt);
time_max = cal_time / dt;
cycle = 1/freq;
time = dt : dt : cal_time;
t = 1 : time_max_g;
pi = 3.1416;
p_keisoku_taihi = zeros(1, time_max_g);
p_keisoku_in = zeros(1, time_max_g);

% f1 = figure;
% figure(f1);
i_x = 1 : ix + 1;
x = i_x * dx;

pin = zeros(size(time));
j = sqrt(-1);

for i = 1 : time_max
    if sweep == 1
        omega = 2*pi*(freq_start + (freq_add)*time(i)/cal_time);
    else
        omega = 2 * pi * freq;
    end
    
    pin(i) = sin(omega * time(i));
end

for t = 1 : time_max_g
    
    for i = 2 : ix
        p1(i) = p2(i) - (rou*c^2*dt/dx)*(u2(i) - u2(i - 1));
    end
    
    p1(1) = 0;
    p2(ix + 1) = 0;
    
    p1(id) = pin(t);
    
    for i = 2 : ix
        p2(i) = p1(i);
    end
    
    for i = 2 : ix
        u1(i) = u2(i) - (dt/(rou*dx))*(p2(i + 1) - p2(i));
    end
    
    p1(end_point) = 0;
    p2(end_point) = 0;
    u2(end_point) = 0;
    u1(end_point) = 0;
    
    if hekomi_bool == 1
        p1(hekomi_point) = 0;
        p2(hekomi_point) = 0;
        u1(hekomi_point) = 0;
        u2(hekomi_point) = 0;
    end
    
    for i = 2 : ix
        u2(i) = u1(i);
    end
    
    p_keisoku_taihi(t) = p1(id + 1);
    p_keisoku_in(t) = p_keisoku_taihi(t) - pin(t);
    
    if mod(t, 50) == 0 
        time(t)
    end
    
    if mode_plot == 0
        if mod(t, 50) == 0
            plot(x(id : end_i), real(p1(id : end_i)));
            grid on;
            drawnow
        end
    end

    if mode_plot == 1
        if mod(t, 50) == 0
            plot(time(1:t), p_keisoku_in(1:t))
            grid on;
            drawnow;
        end
    end
    if t == time_max_g
        dlmwrite('1d1000.csv', p_keisoku_taihi, 'precision', '%.10f', 'delimiter', ',')
        disp("‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ‚ ")
    end
end
