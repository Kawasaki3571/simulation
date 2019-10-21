Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
% M = csvread('二万点データ.csv',2,5);
M = csvread('myFile6000col.csv');
% i = 1 : 17898;
i = 1 : 27000;
data = M(i,1);
plot(i,data)