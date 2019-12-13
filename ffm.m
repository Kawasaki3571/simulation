Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
% M = csvread('二万点データ.csv',2,5);
M = csvread('test.csv');
data = M;
data_size = length(data);
% data_size = size(data);
% i = 1 : 17898;
i = 1 : data_size;
time = i * 0.18 / data_size;
freq = 1000 + i * 4000/data_size;
data_start = 1000;
data_end = 50000;
% plot(freq(500:data_size/4),data(500:data_size/4))
plot(freq(data_start:data_end),data(data_start:data_end));
xlabel("freqency(Hz)")
% plot(i,data)

f2 = figure;
figure(f2)

M2 = csvread('kairyouhanheko0300.csv');
data2 = M2;
data_size2 = length(data2);
% data_size = size(data);
% i = 1 : 17898;
i = 1 : data_size2;
time = i * 0.18 / data_size2;
freq = 1000 + i * 4000/data_size2;
% plot(freq(500:data_size/4),data(500:data_size/4))
plot(freq(data_start:data_end),data2(data_start:data_end));
%plot(freq,data2)
xlabel("freqencyyyyyy(Hz)")
% plot(i,data)