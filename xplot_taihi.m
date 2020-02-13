mode = 0 ;
if mode == 0
    y_10 = csvread('x_0.5_10deg.csv');
    y_20 = csvread('x_0.5_20deg.csv');
    y_30 = csvread('x_0.5_30deg.csv');
    y_40 = csvread('x_0.5_40deg.csv');
    y_50 = csvread('x_0.5_50deg.csv');
    y_60 = csvread('x_0.5_60deg.csv');
    hold on
    plot(x*10^3, y_10)
    plot(x*10^3, y_20)
    plot(x*10^3, y_30)
    plot(x*10^3, y_40)
    plot(x*10^3, y_50)
    plot(x*10^3, y_60)
    legend('10', '20', '30', '40', '50', '60')
    xlim([0 130]);
    hold off
end
if mode== 1
end