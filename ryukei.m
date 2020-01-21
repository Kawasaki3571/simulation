p0 = 1;
r = 1;
cal_time = 0.18;
dt = cal_time / 90000;
t = 0 : dt : cal_time - dt;
freq_start = 1000;
freq_add = 3000;
c = 340;
pai = 3.1415;
sweep = 1;

if sweep == 1
    freq = freq_start + freq_add*t/0.18;
else
    freq_start = 4000;
    freq_add = 0;
    freq = freq_start + freq_add*t/0.18;
end
omega = freq * 2 * pai;
hasu = 2 * pai * freq / c;
j = sqrt(-1);
xL = 0.5;
x = 0 : 0.001 : 2;

p = p0.*(1 + exp(j .* 2 .* hasu .* xL)) .* exp(-1 * j .* omega .* t - pai/2);


% px = p0.*(exp(j .* hasu .* x)) .* exp(-1 * j .* omega .* t - pai/2);

% pt = p0.*(exp(j .* hasu .* x)) .* exp(-1 * j .* omega .* t - pai/2  - (0.6/340));
pr = zeros(1, 90000);
pi = zeros(1, 90000);
po = zeros(1, 90000);

x = 2*xL;
for i = 1 : 90000
    if i * dt >= 2 * xL/340
%         pr(i) = p0.*(exp(j .* hasu(i) .* x)) .* exp(-1 * j .* omega(i) .* t(i) - pai/2  - (2 * xL/340));
        pr(i) = p0.*(exp(j .* hasu(i) * 2.* xL)) .* exp(-1* j .* omega(i) .* t(i)- j*(2 * xL/340));
    else
        pr(i) = 0;
    end
end

x = 0;
for i = 1 : 90000
    if i * dt >= 0
        pi(i) = p0.*(exp(j .* hasu(i) .* x)) .* exp(-1 * j .* omega(i) .* t(i));
    else
        pi(i) = 0;
    end
end

pr  = pr';
pi = pi';
po = pi + pr;

%px = p0.*(exp(j .* hasu .* x)) ;
%p = p0.*(1 + exp(j .* 2 .* hasu .* xL)) ;
p = p';

f2 = figure;
figure(f2);

kasou_out5 = csvread('1d0500kasou5.csv');
kasou_out1 = csvread('1d0500nama.csv');
kasou_out3 = csvread('1d0500kasou3.csv');

kasou_in = csvread('1dnoloadkasou.csv');
nosweep = csvread('1dnosweep.csv');

gousei = kasou_out + pi;
% gousei = kasou_out + pi;
kaishi = 0.015;
shuryo = kaishi + 0.005;

hold on
plot(t(kaishi /dt + 1 : shuryo/dt), kasou_out1(kaishi /dt + 1 : shuryo/dt));
%plot(t(kaishi /dt + 1 : shuryo/dt), nosweep(kaishi /dt + 1 : shuryo/dt));
plot(t(kaishi /dt + 1 : shuryo/dt), pr(kaishi /dt + 1 : shuryo/dt));
hold off
%dlmwrite('riron0500.csv', gousei, 'precision', '%.10f', 'delimiter', ',')
%dlmwrite('rironnoload.csv', pi, 'precision', '%.10f', 'delimiter', ',')
%shuhasu = 1000 + 3000*(kaishi/0.18);
%1/shuhasu
