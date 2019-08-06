ganma = 72.8;
fs = 100000 : 5000000 ;
% fs = 2000000;
pai = 3.1416;
rou = 0.997;
ryukei_d = 0.34 * ((8 * pai * ganma)./(rou * (fs).^2)) .^(1/3);
plot(fs, ryukei_d)
xlabel("fs(Hz)")
ylabel("D(cm)")
title("Relationsip of frequency and diameter of particle")