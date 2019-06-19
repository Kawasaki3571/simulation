xrange = 2.02;
yrange = 0.02;
dx = 0.001;
% dx、発散しないために
% dt = 0.000008;
dt =   0.000001;
% クー数から条件を立てる
cal_time = 0.01;
x1 = 0.03;
x2 = 0.04;
y1 = 0.1;
y2 = 0.44;
xd = 0.008;
xd2 = 0.012;
yd = 0.001;
yd2 = 0.002;
disp_hensu = 0;
f = 1000; %周波数
c0 = 340; %音速
rou0 = 1.293; %密度（kg/m^3
absp0 = - 0.5; % 吸収係数
gensui0 = (f*absp0) / (8.686*c0); % 減衰係数
ix = round(xrange / dx); %x空間感覚の数
jx = round(yrange / dx); %y空間感覚の数
tx = round(cal_time / dt ); %時間感覚の数
%  １おわ
td = 15; % 周波数の代入（？）
id = round(xd / dx);%xの位置（？）
id2 = round(xd2 / dx ) ;
jd = round(yd / dx) ; %y空間感覚の代入（？）
jd2 = round(yd2/dx);
snap = 10 ;%音圧ファイルの出力感覚
del_T = round((cal_time/dt)/snap);%スナップショットの時間感覚の代入
ix1 = round(x1/dx);
ix2 = round(x2/dx);
jx1 = round(y1/dx);
jx2 = round(y2/dx);
p1 = zeros(ix+1,jx+1);
p2 = zeros(ix+1,jx+1);
p3 = zeros(ix+1,jx+1);
u1 = zeros(ix+1,jx+1);
u2 = zeros(ix+1,jx+1);
v1 = zeros(ix+1,jx+1);
v2 = zeros(ix+1,jx+1);
mo = zeros(ix+1,jx+1); %チューブモデル
pressure = zeros(ix + 1 , jx + 1);
SC = 0 ;%励振関数 ０なら連続1ならガウシアン2ハニング3正弦波数波
WN = 1;
W_end = round(2*WN/(f*dt)) - 1 ;
% pin = zeros((round(2*WN/(f*dt))) - 1)
pin = zeros(1,tx)
% crn =(c0 * dt)/dx ; %クーラン数
crn = 0.2;
dd = 199/200 ;%higdons absorption boundary
pai = 3.1415 ;
hasu_o0 = 2*pai*f/c0 ;%実際の波数
hasu0 = sqrt(hasu_o0*hasu_o0 - gensui0^2);%損失ある場合の波数
c_m0 = 2*pai*f/hasu0;%損失のある場合の音速
alpha0 = 2*hasu_o0*rou0*c_m0/hasu0;%吸収項
kap0 = c_m0^2*rou0;
cp1 = 1;
cp2 = kap0*dt/dx;
cv1 = (2*rou0 - alpha0*dt)/(2*rou0 + alpha0*dt);
cv2 = dt*2 /((2*rou0 + alpha0*dt)*dx);
a = 2*(crn - 1)/(crn + 1);
b = ((crn - 1)/(crn + 1))*2;%要検討要検討
c = 2 * (crn - 1) / (crn + 1) * dd;
d = dd * 2;
e = dd *2; %要検討要検討
a1 = 1 / cos(0);
a2 = 1 / cos(0);
d1 = 0.005;
d2 = 0.005;
ca0 = (c_m0 * dt - dx) / ((c_m0 * dt + dx));
ca1 = (dx * 2) / (c_m0 * dt + dx);
ca2 = (dx * c_m0 * dt * c_m0 * dt) / (dx * dx * 2 * (c_m0 * dt + dx));
cah1 = (a1 * c_m0 * dt - dx) / (a1 * c_m0 * dt + dx);
cah2 = (a2 * c_m0 * dt - dx) / (a2 * c_m0 * dt + dx);
a = cah1 + cah2;
b = cah1 * cah2;
c = cah1 * (1 - d2) + cah2 * (1 - d1);
d = ((1 - d1) + (1 - d2));
e = ((1 - d1) * (1 - d2));
%%初期化①
x=0:dx:xrange; 
y=0:dx:yrange;
for i = jx1 : jx2
    for j = jx1 : jx2
            mo(i,j) = 0 ;
    end
end
sn = 0;
w = 2 * pai * f ;
switch SC
    case 0
        for t = 1 : W_end
            tr = w * dt * t / WN;
            if tr <= pai
                pin(t) = sin((w * dt * t));
            else
                pin(t) = sin(w * dt * t);
            end
        end
    case 1
        tau = WN / f;
        for t = 1:W_end
            al = (4/tau) * (4/tau);
            pin =exp(-al * (dt*t - tau)*(dt*t - tau)) * sin(w*(dt*t - tau));
        end
            %%②
    case 2
        for t = 1 : W_end
            pin = (1-cos(w*dt*t/WN)) * sin(w*dt*t)/2;
        end

        %%②
    case 3
        for t = 1 : W_end
        pin = sin(w * dt * real(t - 1));
        end
end
%%timeloop
for t = 1: tx
   time = t * dt;
    for i = 2:ix
        for j = 2:jx
            p1(i,j) = cp1 .* p2(i,j) - cp2 .* (u2(i,j) - u2(i-1,j) + v2(i,j) - v2 (i,j-1));
        end
    end



    i = 1;
    for j = 2 : jx
        p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
    end
    i = ix + 1;
    for j = 2 : jx
        p1(i,j) = a .* (p1(i - 1,j) - p2(i,j)) - b .* (p1(i-2,j) - 2*p2(i-1,j) + p3(i,j)) - c .* (p2(i - 2,j) - p3(i-1,j)) + d .* p2(i -1,j) - e .* p3(i-2,j);
    end
%     ⑤おわ
    j = 1;
    for i = 2 : ix
        p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
    end
    j = jx + 1;
    for i = 2 : ix
        p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
    end
%     if t <= tx
         for i = id : id2
            for j = jd : jd2
                p1(i,j) = pin(t);
            end
         end
%      end
% ああああああ
    %soundpressureのスワップ
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            p3(i,j) = p2(i,j);
            p2(i,j) = p1(i,j);
        end
    end
    %xのvelocityの更新
    for i = 1 : ix
        for j = 2 : jx
            u1(i,j) = cv1 .* u2(i,j) - cv2 .* (p2(i + 1, j) - p2(i,j));
        end
    end
    %y
    for i = 2 : ix
        for j = 1 : jx
            v1(i,j) = cv1 * v2(i,j) - cv2 * (p2(i, j + 1) - p2(i,j));
        end
    end
    for i = 1 : ix
        for j = 1 : jx
            if mo(i,j) ~= mo(i + 1,j) then
                u1(i,j) = 0;
            end
            if mo(i,j) ~= mo(i,j + 1) then
                v1(i,j) = 0;
            end
        end
    end
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            u2(i,j) = u1(i,j);
            v2(i,j) = v1(i,j);
        end
    end
    %④
    if  mod(t, snap) == 0
        sn = sn + 1;
    for j = 1 : jx
        for i = 1 ; ix
%             p2(i,j);
        end
    end
    end
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            pressure(i,j) = (p1(i,j) .^2).^0.5 ;
        end
    end
%        if p1(200,1) > 0
%            disp("0.2地点がプラスです！！！");
%        end
%        if p1(400,1) > 0
%            disp("0.4地点がプラスです");
%        end
%     if p1(200, 10) > 0
%         disp_hensu = 1;
%     end
%     if p1(400,10) > 0
%         disp_hensu = 2;
%     end
%     if disp_hensu == 1
%         disp("0.2地点がプラスになりました！！")
%         time
%     end
%     if disp_hensu == 2
%         disp("0.4地点がプラスになりました！！")
%         time
%     end
    y = 1 : jx + 1 ;
    x = 1 : ix + 1;
%    [X,Y] = meshgrid(x,y);
%    Z = zeros(ix + 1, jx + 1);
%    for i = 1:ix + 1
%       for j = 1:jx + 1
%           Z(i,j) = pressure(X,Y);
%       end
%    end
%    mesh(X,Y,Z);
%    surf(X,Y,Z)
       imagesc(y,x,pressure(x,y))
       title(['pressure when ',num2str(time),'seconds have passed'])
       xlabel('x(mm)')
       ylabel('y(mm)')
%    mesh(y,x,pressure);
    grid on;
    drawnow
end
% jjjaa