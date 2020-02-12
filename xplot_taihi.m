y_20 = csvread('x_20deg.csv');
y_40 = csvread('x_40deg.csv');
hold on
plot(x*10^3, y_20)
plot(x*10^3, y_40)
legend('20', '40')
xlim([0 130]);
hold off