p0 = 1;
r = 1;
cal_time = 0.18;
dt = cal_time / 90000;
t = 0 : dt : cal_time;
freq_start = 1000;
freq_add = 3000;
c = 340;
pai = 3.1415;

freq = freq_start + freq_add*t/0.18;
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

% plot(t(1 : 0.18/dt), po(1 : 0.18/dt)); 
hold on
% plot(t(1 : 0.02/dt), pr(1 : 0.02/dt));
% plot(t(1 : 0.02/dt), pi(1 : 0.02/dt)); 
plot(t(1 : 0.02/dt), po(1 : 0.02/dt)); 
hold off
% dlmwrite('riron0500.csv', po, 'precision', '%.10f', 'delimiter', ',')
% dlmwrite('rironnoload.csv', pi, 'precision', '%.10f', 'delimiter', ',')