Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
% M = csvread('二万点データ.csv',2,5);
M = csvread('二万点データ.csv',2,1);
% i = 1 : 17898;
i = 1 : 17888;
time = M(i,1);
diff = M(i,4);
Y = fft(diff);
P2 = abs(Y/L);
disp(diff)