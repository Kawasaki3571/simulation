mode_plot = 10;
hekomi_bool = 0;
sweep = 1;
input_mode = 0;

c = 340;
rou = 1.2;
baf = 2.00;
l = 2 + baf;
x_in = 0;
hekomi = 0.3 + baf;
keisokuten = x_in + baf;
freq = 2000;
freq_start = 1000;
freq_add = 3000;
ramuda = c / freq ;
dx = ramuda / 100;
crn = 0.4;

dt = dx * crn / c;
ix = round(l / dx);
id = round(x_in / dx) + 1;
sokutei_end_x  = l;
end_i = sokutei_end_x / dx;
end_point = round(end_i);
hekomi_point = round(hekomi/dx);
keisokuten = round(keisokuten/ dx);
p1 = zeros(1, ix + 1);
p2 = zeros(1, ix + 1);
p3 = zeros(1, ix + 1);
u1 = zeros(1, ix + 1);
u2 = zeros(1, ix + 1);
cal_time = 0.18;

dd_b = 199/200;
ddd = 999/1000;
%ddd = 5001/5000;
%ddd = 1;
a_b = 2 * (crn-1)/(crn+1);
b_b = ((crn - 1) / (crn + 1)) ^ 2;
c_b = ((crn - 1) / (crn + 1)) * dd_b*2;
d_b = dd_b * 2;
e_b = dd_b ^ 2;

time_max_g = round(cal_time/dt);
time_max = cal_time / dt;
cycle = 1/freq;
time = dt : dt : cal_time;
t = 1 : time_max_g;
pi = 3.1416;
p_keisoku_taihi = zeros(1, time_max_g);
p_keisoku_taihi2 = zeros(1, time_max_g);
p_keisoku_spec = zeros(1, time_max_g);
p_keisoku_in = zeros(1, time_max_g);
p_keisoku_shin = zeros(1, time_max_g);

% f1 = figure;
% figure(f1);
i_x = 1 : ix + 1;
x = i_x * dx;

pin = zeros(size(time));
pin_gyaku = zeros(size(time));
j = sqrt(-1);

if sweep == 1
    %frequency = freq_start + (freq_add)*time(i)/cal_time;
else
    frequency = freq;
end

pin_range = time_max / 1;
if input_mode == 0
    for i = 1 : pin_range
        if sweep == 1
            omega = 2*pi*(freq_start + (freq_add)*time(i)/cal_time);
        else
            omega = 2 * pi * freq;
        end
        pin(i) = 3 * exp(-1 * j * omega * time(i) - (j * pi/2));
    end
end

if input_mode == 1
    for i = 1 : 10
        pin(i) = 1;
    end
end

if input_mode == 2
    if sweep == 0
        cycle_point = round(cycle/ dt);
        omega = 2 * pi * freq;
        for i = 1 : cycle_point
            pin(i) = exp(j * omega * time(i) - (j * pi/2));
        end
    else
        pin = 0;
    end
end

for i = 1 : time_max
    pin_gyaku(i) = 1 ./ pin(i);
end

for t = 1 : time_max_g

    for i = 2 : ix
        p1(i) = ddd*(p2(i) - (rou*c^2*dt/dx)*(u2(i) - u2(i - 1)));
    end


    %i = 1;
    i = 1;
    p1(i) = a_b .* (p1(i + 1)-p2(i)) - b_b .* (p1(i + 2) - 2 .* p2(i + 1) + p3(i)) - c_b .* (p2(i + 2) - p3(i + 1)) + d_b .* p2(i + 1) - e_b .* p3(i + 2);
    i = end_point;
    p1(i) = a_b .* (p1(i - 1)-p2(i)) - b_b .* (p1(i - 2) - 2 .* p2(i - 1) + p3(i)) - c_b .* (p2(i - 2) - p3(i - 1)) + d_b .* p2(i - 1) - e_b .* p3(i - 2);
    
    if t < pin_range
        p1(id) = pin(t);
    end
    
    for i = 1 : ix
        p3(i) = p2(i);
        p2(i) = p1(i);
    end
    
    for i = 1 : ix
        u1(i) = u2(i) - (dt/(rou*dx))*(p2(i + 1) - p2(i));
    end
    
    %p1(end_point) = 0;
    %p2(end_point) = 0;
    u2(end_point) = 0;
    u1(end_point) = 0;
    if hekomi_bool == 1
        u1(hekomi_point) = 0;
        u2(hekomi_point) = 0;
    end
    
    for i = 1 : ix
        u2(i) = u1(i);
    end
    
    p_keisoku_taihi(t) = p1(keisokuten + 1);
    p_keisoku_taihi2(t) = p1(round(hekomi_point / 2));
    p_keisoku_in(t) = p_keisoku_taihi(t) - pin(t);
    p_keisoku_spec(t) = abs(p_keisoku_taihi(t)).^2;
    
    p_keisoku_shin(t) = p_keisoku_taihi(t) .* pin_gyaku(t);
    
    if mod(t, 50) == 0 
        time(t)
    end
    
    if mode_plot == 0
        if mod(t, 50) == 0
            plot(x(keisokuten : round(l/dx)), real(p1(keisokuten : round(l/dx))));
            grid on;
            drawnow
        end
    end

    if mode_plot == 1
        if mod(t, 50) == 0
%             plot(time(1:t), p_keisoku_spec(1:t))
            plot(time(1:t) - 2/340, p_keisoku_taihi(1:t))
            grid on;
            drawnow;
        end
    end
    
    if mode_plot == 2
        if mod(t, 50) == 0
            plot(time(1:t), p_keisoku_shin(1:t))
            grid on;
            drawnow;
        end
    end
    
    if mode_plot == 3 && hekomi_bool == 0
        if mod(t, 50) == 0
            plot(time(1:t), p_keisoku_taihi2(1:t))
            grid on;
            drawnow;
        end
    end
    
    if t == time_max_g
        p_keisoku_taihi = p_keisoku_taihi';
        p_keisoku_shin = p_keisoku_shin';
        %p_keisoku_spec = p_keisoku_spec';
        dlmwrite('1d1000.csv', p_keisoku_taihi, 'precision', '%.10f', 'delimiter', ',')
        disp("あああああああああああああああああ")
    end
end
