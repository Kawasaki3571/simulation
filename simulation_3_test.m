xrange = 0.02;
yrange = 0.02;
zrange = 0.02;
dx = 0.001;
dt =   0.000002;

cal_time = 0.01;

mode = 0; %0...通常シミュレーション 1...初期印加位置表示
mode_plot = 0; %プロットモード選択 0...カラーマップ進行 1...xプロット 2...xプロット進行

yd = 0.003;
yd2 = 0.017;
zd = 0.003;
zd2 = 0.017;
xd = 0.001;
xd2 = 0.002;
t1 = 0;
t2 = 0;
speed = 0;
disp_hensu = 0;
f = 1000; %周波数
c0 = 340; %音速
rou0 = 1.293; %密度（kg/m^3
absp0 = - 0.5; % 吸収係数
gensui0 = (f*absp0) / (8.686*c0); % 減衰係数
ix = round(xrange / dx); %x空間感覚の数
jx = round(yrange / dx); %y空間感覚の数
kx = round(zrange / dx) ;
tx = fix(cal_time / dt ); %時間感覚の数
%  １おわ
td = 15; % 周波数の代入（？）
id = round(xd / dx);%xの位置（？）
id2 = round(xd2 / dx ) ;
jd = round(yd / dx) ; %y空間感覚の代入（？）
jd2 = round(yd2/dx);
kd = round(zd / dx);
kd2 = round(zd2 / dx);
snap = 10 ;%音圧ファイルの出力感覚
del_T = round((cal_time/dt)/snap);%スナップショットの時間感覚の代入
p1 = zeros(ix+1,jx+1,kx+1);
p2 = zeros(ix+1,jx+1,kx+1);
p3 = zeros(ix+1,jx+1,kx+1);
u1 = zeros(ix+1,jx+1,kx+1);
u2 = zeros(ix+1,jx+1,kx+1);
v1 = zeros(ix+1,jx+1,kx+1);
v2 = zeros(ix+1,jx+1,kx+1);
w1 = zeros(ix+1,jx+1,kx+1);
w2 = zeros(ix+1,jx+1,kx+1);
mo = zeros(ix+1,jx+1,kx+1); %チューブモデル
pressure = zeros(ix + 1, jx + 1, kx + 1);
SC = 0 ;%励振関数 ０なら連続1ならガウシアン2ハニング3正弦波数波
WN = 1;
W_end = round(2*WN/(f*dt)) - 1 ;
pin = zeros(1,tx);
p_keisoku_taihi = zeros(1,tx);
crn =(c0 * dt)/dx ; %クーラン数
dd = 199/200 ;%higdons absorption boundary
pai = 3.1415 ;
hasu_o0 = 2*pai*f/c0 ;%実際の波数
hasu0 = sqrt(hasu_o0*hasu_o0 - gensui0^2);%損失ある場合の波数
c_m0 = 2*pai*f/hasu0;%損失のある場合の音速
alpha0 = 2*hasu_o0*rou0*c_m0/hasu0;%吸収項
kap0 = c_m0^2*rou0;
cp1 = 1;
cp2 = rou0*c_m0^2*dt/dx;
cv1 = 1;
cv2 = dt / (rou0 * dx)
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
for i = 1 : ix + 1
    for j = 1 : jx + 1
        for k = 1 : kx + 1
            p1(i,j,k) = 0;
            p2(i,j,k) = 0;
            p3(i,j,k) = 0;
            u1(i,j,k) = 0;
            u2(i,j,k) = 0;
            v1(i,j,k) = 0;
            v2(i,j,k) = 0;
            p1_taihi(i,j) = 0;
        end
    end
end
x=0:dx:xrange; 
y=0:dx:yrange;
z=0:dx:zrange;

sn = 0;
w = 2 * pai * f ;
switch SC
    case 0
        for t = 1 : tx 
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
    case 2
        for t = 1 : W_end
            pin = (1-cos(w*dt*t/WN)) * sin(w*dt*t)/2;
        end
    case 3
        for t = 1 : W_end
        pin = sin(w * dt * real(t - 1));
        end
end

if xd > xd2
    error("xrangeがおかしい")
end
if yd > yd2
    error("yrangeがおかしい")
end
if zd > zd2
    error("zrangeがおかしい")
end
if id <= dx || id2 > ix || jd <= dx || jd2 > jx
    error("rangeおかしい");
end

if mode == 0
%timeloop
for t = 1: tx
   time = t * dt * 1.47;
   if crn >= 1
       error("クーラン数が不適切です。");
       break;
   end
    for i = 2:ix
        for j = 2:jx
            for k = 2:kx
                if t <= tx && i >= id && i <= id2 && j >= jd && j <= jd2 && k >= kd && k <= kd2
                    p1(i,j,k) = pin(t);
                    p2(i,j,k) = pin(t);
                else
                    p1(i,j,k) = cp1 .* p2(i,j,k) - cp2 .* (u2(i,j,k) - u2(i-1,j,k) + v2(i,j,k) - v2 (i,j-1,k) + w2(i,j,k) - w2(i,j,k-1));
                end
            end
        end
    end
    for k = 1 : kx + 1
    i = 1;
    for j = 2 : jx
        p1(i,j,k) = a .* (p1(i + 1,j,k)-p2(i,j,k)) - b .* (p1(i + 2,j,k) - 2 .* p2(i + 1,j,k) + p3(i,j,k)) - c .* (p2(i + 2,j,k) - p3(i + 1,j,k)) + d .* p2(i + 1,j,k) - e .* p3(i + 2,j,k);
%         p1(i + 1, j) = p1_taihi(i + 1, j);
    end
    i = ix + 1;
    for j = 2 : jx
        p1(i,j,k) = a .* (p1(i - 1,j,k) - p2(i,j,k)) - b .* (p1(i-2,j,k) - 2*p2(i-1,j,k) + p3(i,j,k)) - c .* (p2(i - 2,j,k) - p3(i-1,j,k)) + d .* p2(i -1,j,k) - e .* p3(i-2,j,k);
         p1(i - 1, j) = p1_taihi(i - 1, j);
    end
    j = 1;
    for i = 2 : ix
        p1(i,j,k) = a .* (p1(i,j+1,k) - p2(i,j,k)) - b .* (p1(i,j+2,k) - 2 * p2(i,j+1,k) + p3(i,j,k)) - c .* (p2(i,j+2,k) - p3(i,j + 1,k)) + d .* p2(i,j + 1,k) - e .* p3(i,j + 2,k);
%         p1(i ,j + 1) = p1_taihi(i, j + 1);
    end
    j = jx + 1;
    for i = 2 : ix
        p1(i,j,k) = a .* (p1(i,j-1,k) - p2(i,j,k)) - b .* (p1(i,j-2,k) - 2 * p2(i,j-1,k) + p3(i,j,k)) - c .* (p2(i,j-2,k) - p3(i,j - 1,k)) + d .* p2(i,j - 1,k) - e .* p3(i,j - 2,k);
%         p1(i, j - 1) = p1_taihi(i,j - 1);
    end
    end
    
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            for k = 1 : kx + 1
                p3(i,j,k) = p2(i,j,k);
                p2(i,j,k) = p1(i,j,k);
            end
        end
    end
    for i = 1 : ix
        for j = 2 : jx
            for k = 1 : kx + 1
                u1(i,j,k) = cv1 .* u2(i,j,k) - cv2 .* (p2(i + 1, j,k) - p2(i,j,k));
            end
        end
    end
    for i = 2 : ix
        for j = 1 : jx
            for k = 1 : kx + 1
                v1(i,j,k) = cv1 * v2(i,j,k) - cv2 * (p2(i, j + 1,k) - p2(i,j,k));
            end
        end
    end
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            for k = 1 : kx + 1
                u2(i,j) = u1(i,j);
                v2(i,j) = v1(i,j);
            end
        end
    end
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            for k = 1 : kx + 1
                pressure(i,j,k) = (p1(i,j,k) .^2).^0.5 ;
            end
        end
    end
%     if t1 == 0
%         if p1(200,11) > 0
%             t1 = time;
%         end
%     end
%     if t2 == 0
%         if p1(400,11) > 0
%             t2 = time;
%         end
%     end
    
    p_keisoku_taihi(t) = p1(10,10);
    disp(t);
    y = 1 : jx + 1 ;
    x = 1 : ix + 1;
    z = 1 : kx + 1;
    if mode_plot == 0
        if mod(t,10) == 0
            scatter3(x,y,z,pressure(x,y,z));
%         imagesc(y,x,z,pressure(x,y,z));
%           colorbar
            title(['pressure when ',num2str(time),'seconds have passed'])
            xlabel('y(mm)')
            ylabel('x(mm)')
            zlabel('z(mm)')
        grid on;
        drawnow
        end
    end
    if mode_plot == 1
        if t == 5000
            plot(x,p1(x, 10));
            break;
        end
    end
    if mode_plot ==2
        if mod(t,10) == 0
        plot(x,pressure(x, 10));
        grid on;
        drawnow
        end
    end
    if mode_plot == 3
        if t == 5000
            t_x = 1 : t;
            plot(t_x,p_keisoku_taihi(t_x));
            break;
        end
    end
end
end
if mode == 1
    for i = id : id2
        for j = jd : jd2
            p1(i,j) = 100;
        end
    end
     for i = 1 : ix + 1
         for j = 1 : jx + 1
             pressure(i,j) = (p1(i,j) .^2).^0.5 ;
         end
     end
    
    x = 1 : ix + 1 ;
    y = 1 : jx + 1;
    imagesc(y,x,pressure(x,y));
    colorbar;
       title(['Mode 1'])
       xlabel('y(mm)')
       ylabel('x(mm)')
    grid on;
    drawnow
    disp("クーラン数は" + crn );
    disp(size(p1));
    disp(size(pressure));
end