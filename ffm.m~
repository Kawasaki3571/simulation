Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

X = S + 2*randn(size(t));

plot(1000*t(1:L),X(1:L))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

Y = fft(X);

P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
P1 = P2;
P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
f = Fs*(0:L-1)/L;
% disp(size(P1))
% disp(size(f))

plot(f,(P1))
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


