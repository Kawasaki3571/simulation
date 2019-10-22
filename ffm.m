Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
data_size = 17900;
% M = csvread('二万点データ.csv',2,5);
M = csvread('test.csv');
% i = 1 : 17898;
i = 1 : data_size;
time = i * 0.18 / 17900;
freq = 2000 + i * 7000/data_size;
data = M(i,1);
plot(freq,data)
% plot(i,data)