
% syms a b c x
% eqn = a*x^2 + b*x + c == 0;
% solx = solve(eqn, x)
theta = pi:0.01:3*pi;
n = 100;
d_d = 0.002;
b_haba = d_d / n;
syms x
dhi1 = (sin(n*3.14*b_haba*freq*sin(theta + 0.5*pi)/340)) ./(n * sin(3.14*b_haba*freq*sin(theta + 0.5*pi)/340)) == 0.5 ;
solx = solve(dhi1, theta);
%     if mod(freq,5000) == 0
%         polarplot(theta,abs(dhi1));
%         title(['Directivity when frequency is',num2str(freq),'Hz'])
%         drawnow
%         pause(1);
%     end