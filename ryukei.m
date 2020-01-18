p0 = 1;
r = 1;
cal_time = 0.18;
t = 0 : dt : cal_time;
freq_start = 1000;
freq_add = 3000;
c = 340;
freq = freq_start + freq_add*t/0.18;
hasu = 2 * pi * freq / c;
j = sqrt(-1);
xL = 0.8;
p = 1 + exp(j .* 2 .* hasu .* xL);
plot(t , abs(p))