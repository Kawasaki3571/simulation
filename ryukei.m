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
gosahyoka = 0;
csvrangemax = 90000;

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
xL = 1.0;
x = 0 : 0.001 : 2;

p = p0.*(1 + exp(j .* 2 .* hasu .* xL)) .* exp(-1 * j .* omega .* t - pai/2);


% px = p0.*(exp(j .* hasu .* x)) .* exp(-1 * j .* omega .* t - pai/2);

% pt = p0.*(exp(j .* hasu .* x)) .* exp(-1 * j .* omega .* t - pai/2  - (0.6/340));
pr = zeros(1, 90000);
pr0 = zeros(1, 90000);
pr2 = zeros(1, 90000);
pi = zeros(1, 90000);
po = zeros(1, 90000);

x = 2*xL;
counter = 0;
for i = 1 : 90000
    if i * dt >= 2 * xL/340
%         pr(i) = p0.*(exp(j .* hasu(i) .* x)) .* exp(-1 * j .* omega(i) .* t(i) - pai/2  - (2 * xL/340));
        if counter == 0
            i_hasu = i-1;
        end
        counter = 1;
        pr(i) = p0.*(exp(j .* hasu(i) * 2.* xL)) .* exp(-1* j .* omega(i) .* t(i)- j*(2 * xL/340));
        %pr2(i) = p0.*(exp(j .* hasu(i - i_hasu) * 2.* xL)) .* exp(-1* j .* omega(i - i_hasu) .* t(i)- j*(2 * xL/340));
        pr2(i) = p0.*(exp(j .* hasu(i - i_hasu) * 2.* xL)) .* exp(-1* j .* omega(i - i_hasu) .* t(i)- j*(2 * xL/340));
    else
        pr(i) = 0;
    end
end

for i = 1 : 90000
    %pr2(i) = p0.*(exp(j .* hasu(i) * 2.* xL)) .* exp(-1* j .* omega(i) .* t(i));
    pr0(i) = p0.*(exp(j .* hasu(i) * 2.* xL)) .* exp(-1* j .* omega(i) .* t(i));
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
pr2 = pr2';
pr0 = pr0' + 0.01;
pi = pi';
po = pi + pr;
po2 = pi + pr2;

%px = p0.*(exp(j .* hasu .* x)) ;
%p = p0.*(1 + exp(j .* 2 .* hasu .* xL)) ;
p = p';

f2 = figure;
figure(f2);

kasou_out1 = csvread('1d0500nama.csv');
kasou_out2 = csvread('1d0500kasou.csv');
kasou_out3 = zeros(90000, 1);

achieve_time = (2*xL)/c0;
achieve_time_g = round(achieve_time / dt);

kasou_out3(1 : csvrangemax - achieve_time_g + 1) = kasou_out1(achieve_time_g : csvrangemax);

kasou_in = csvread('1dnoloadkasou.csv');
nosweep = csvread('1dnosweep.csv');

gousei = kasou_out3 + pi;
% gousei = kasou_out + pi;
kaishi = 0.0;
shuryo = kaishi + 0.18;

kasou_outnyu = zeros(1, round((cal_time-2*xL/340)/dt));
kasou_outnyu = kasou_out2(round((2*xL/340)/dt) : round(cal_time/dt));
pr_outnyu = pr(round((2*xL/340)/dt) : round(cal_time/dt));

hold on
%plot(t(kaishi /dt + 1 : shuryo/dt), kasou_out3 (kaishi /dt + 1 : shuryo/dt));
plot(t(kaishi /dt + 1 : shuryo/dt), po2(kaishi /dt + 1 : shuryo/dt));
%plot(t(kaishi /dt + 1 : shuryo/dt), kasou_out1(kaishi /dt + 1 : shuryo/dt));
%legend('1éü', '2éü', 'óùò_íl')
hold off

% plot(t, po)
dlmwrite('riron1000zure.csv', po2, 'precision', '%.10f', 'delimiter', ',')
dlmwrite('rironnoloadzure.csv', pi, 'precision', '%.10f', 'delimiter', ',')

