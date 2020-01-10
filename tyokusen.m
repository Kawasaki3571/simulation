c = 340;
rou = 1.2;
l = 100;
x_in = 50;
freq = 100;
ramuda = c / freq ;
dx = ramuda / 100;
crn = 0.4;
dt = dx * crn / c;
ix = round(l / dx);
id = round(x_in / dx);
sokutei_end_x  = 55;
end_i = sokutei_end_x / dx;
p1 = zeros(1,ix + 1);
p2 = zeros(1,ix + 1);
u1 = zeros(1, ix + 1);
u2 = zeros(1,ix + 1);
cal_time = 1;
time_max = cal_time / dt;
cycle = 1/freq;
time = dt : dt : cal_time;
t = 1 : round(cal_time / dt);
pi = 3.1416;
omega = 2*pi*freq;

f1 = figure;
figure(f1);

% pin
pin = zeros(size(time));
j = sqrt(-1);
for i = 1 : size(t)
    pin(i) = sin(omega * time(i));
end

for t = 1 : round(cal_time / dt)
    
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
    
    for i = 2 : ix
        u2(i) = u1(i);
    end
    
    i_x = 1 : ix + 1;
    x = i_x * dx;
    
    plot(x(id : end_i), real(p1(id : end_i)));
    
    grid on;
    drawnow
    time(t)
end
