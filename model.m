f1 = figure;
% f2 = figure;

c = 340; %âπë¨
c0 = 340;
% rou0 = 1.293; %ñßìxÅikg/m^3
rou0 = 1.293;
% rou0 = 1000;?
freq_param = 0.8;
freq_a = 3000;
freq_start = 2000;
freq_add = 7000;
% freq = freq_a*freq_param;

freq = freq_a;

freq_abs = freq_a;
ramuda = c0 / freq;
dx_param = 0.05; %0.05-0.025

dx_param = 0.02;

dx = ramuda*dx_param; % É…ÇÃ20-30ï™ÇÃàÍ
% dt = dx / (5 * c);
% dt = dx*0.15 / (c0);
crn_param = 0.2;
dt = dx*crn_param/ (c0);

sei = csvread('kyokaisei.csv');
naname = csvread('kyokainaname.csv');
time = (1 : size(sei)) * dt;
hold on
plot(time(1 : 1000), sei(1 : 1000))
plot(time(1 : 1000), naname(1 : 1000))
hold off