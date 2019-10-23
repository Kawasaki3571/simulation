Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
data_size = 23900;
% M = csvread('二万点データ.csv',2,5);
M = csvread('test.csv');
data = M(i,1);
% data_size = size(data);
% i = 1 : 17898;
i = 1 : data_size;
time = i * 0.18 / data_size;
freq = 2000 + i * 7000/data_size;
plot(freq(500:data_size/4),data(500:data_size/4))
% plot(i,data)