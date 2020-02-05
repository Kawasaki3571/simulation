i = 1 : ix + 1;
rangeA = round(0.1/dx) : round(0.8/dx); 
x_i = i(rangeA) * dx;
px_horn = csvread('horn_x.csv');
px_nothorn = csvread('nothorn_x.csv');
hold on
plot(x_i, px_horn(rangeA))
plot(x_i, px_nothorn(rangeA))
legend('horn', 'tube')
hold off