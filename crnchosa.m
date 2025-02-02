range = 2;
boundary = 2;
mode = 0; %0...通常シミュレーション 1...初期印加位置表示 2...変数表示3...指向性の極グラフ
mode_plot = 1; %プロットモード選択 0...カラーマップ進行 1...xプロット 2...xプロット進行 3...ある地点の時間変化
hekomi = 0;

freq = 5000; %周波数
c = 340;
ramuda = c / freq;
dx = ramuda/25;
% dx、発散しないために
% dt =   0.000002;
dt = dx / (5 * c);
% クー数から条件を立てる
cal_time = 0.01;
if range == 0
    xrange = 0.52;
    yrange = 0.02;
    yd = 0.002;
    yd2 = 0.018;
    xd = 0.002;
    xd2 = 0.003;
end
if range == 1
    xrange = 0.2;
    yrange = 0.2;
    xd = 0.099;
    xd2 = 0.101;
    yd = 0.099;
    yd2 = 0.101;
end
if range == 2
    b_haba = 0.002;
    xrange = 0.5;
    yrange = 0.5;
    xd = 0.249;
    xd2 = 0.251;
    yd_a = 0.25;
    yd_a2 = 0.25;
    yd_b = 0.25 + b_haba;
    yd_b2 = 0.25 + b_haba;
    yd_c = 0.25 + 2*b_haba;
    yd_c2 = 0.25 + 2*b_haba;
    yd_d = 0.25 - b_haba;
    yd_d2 = 0.25 - b_haba;
    yd_e = 0.25 - 2*b_haba;
    yd_e2 = 0.25 - 2*b_haba;
end
if range == 3
    b_haba = 0.002;
    xrange = 0.52;
    yrange = 0.02;
    xd = 0.002;
    xd2 = 0.003;
    yd_a = 0.01;
    yd_a2 = 0.01;
    yd_b = 0.01 + b_haba;
    yd_b2 = 0.01 + b_haba;
    yd_c = 0.01 + 2*b_haba;
    yd_c2 = 0.01 + 2*b_haba;
    yd_d = 0.01 - b_haba;
    yd_d2 = 0.01 - b_haba;
    yd_e = 0.01 - 2*b_haba;
    yd_e2 = 0.01 - 2*b_haba;
end
    
x1 = 0.03;
x2 = 0.04;
y1 = 0.1;
y2 = 0.44;
t1 = 0;
t2 = 0;
speed = 0;
disp_hensu = 0;
c0 = 340; %音速
rou0 = 1.293; %密度（kg/m^3
absp0 = - 0.5; % 吸収係数
b_po = 0.3 ; %凹み位置
h = 0.01;%凹み幅
w = 0.001;%凹みふかさ
gensui0 = (freq*absp0) / (8.686*c0); % 減衰係数
ix = round(xrange / dx); %x空間感覚の数
jx = round(yrange / dx); %y空間感覚の数
tx = fix(cal_time / dt ); %時間感覚の数
b_x = round(b_po / dx);
h_x = round(h / dx);
w_x = round(w / dx);
%  １おわ
td = 15; % 周波数の代入（？）
id = round(xd / dx);%xの位置（？）
id2 = round(xd2 / dx ) ;
if range ~= 2
    jd = round(yd / dx) ; %y空間感覚の代入（？）
    jd_2 = round(yd2/dx);
else
    jd_a = round(yd_a / dx) ; %y空間感覚の代入（？）
    jd_a2 = round(yd_a2/dx);
    jd_b = round(yd_b / dx) ; %y空間感覚の代入（？）
    jd_b2 = round(yd_b2/dx);
    jd_c = round(yd_c / dx) ; %y空間感覚の代入（？）
    jd_c2 = round(yd_c2/dx);
    jd_d = round(yd_d / dx) ; %y空間感覚の代入（？）
    jd_d2 = round(yd_d2/dx);
    jd_e = round(yd_e / dx) ; %y空間感覚の代入（？）
    jd_e2 = round(yd_e2/dx);
end


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
W_end = round(2*WN/(freq*dt)) - 1 ;
pin = zeros(1,tx);
p_keisoku_taihi = zeros(1,tx);
crn =(c0 * dt)/dx ; %クーラン数
dd = 199/200 ;%higdons absorption boundary
pai = 3.1415 ;
hasu_o0 = 2*pai*freq/c0 ;%実際の波数
hasu0 = sqrt(hasu_o0*hasu_o0 - gensui0^2);%損失ある場合の波数
c_m0 = 2*pai*freq/hasu0;%損失のある場合の音速
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
        p1(i,j) = 0;
        p2(i,j) = 0;
        p3(i,j) = 0;
        u1(i,j) = 0;
        u2(i,j) = 0;
        v1(i,j) = 0;
        v2(i,j) = 0;
        p1_taihi(i,j) = 0;
    end
end
x=0:dx:xrange; 
y=0:dx:yrange;
for i = jx1 : jx2
    for j = jx1 : jx2
            mo(i,j) = 0 ;
    end
end
sn = 0;
w = 2 * pai * freq ;
switch SC
    case 0
        for t = 1 : tx 
            time = t * dt;
            tr = w * dt * t;
            if t < pai / (w * dt)
                pin(t) = sin(w * time);
            else
                 pin(t) = sin(w * time);
%                pin(t) = 0;
            end
        end
    case 1
        tau = WN / freq;
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
%%timeloop
if xd > xd2
    error("xrangeがおかしい")
end
% if yd > yd2
%     error("yrangeがおかしい")
% end
if range ~= 2
    if id < 1 || id2 > ix || jd < 1 || jd_2 > jx
    id
    dx
    id2
    ix
    jd
    jd_2
    jx
        error("rangeおかしい");
    end
end

if mode == 0
for t = 1: tx
   time = t * dt;
   if crn >= 1
       disp("クーラン数が不適切です。");
       break;
   end
    for i = 2:ix
        for j = 2:jx
            if range ~= 2
                if t <= tx && i >= id && i <= id2 && j >= jd && j <= jd_2
                    p1(i,j) = pin(t);
                    p2(i,j) = pin(t);
                else
                    p1(i,j) = cp1 .* p2(i,j) - cp2 .* (u2(i,j) - u2(i-1,j) + v2(i,j) - v2 (i,j-1));
                end
            else
                if t <= tx
                    if i >= id && i <= id2
                        if j >= jd_a && j <= jd_a2
                            p1(i,j) = pin(t);
                            %p2(i,j) = pin(t);
                        elseif j >= jd_b && j <= jd_b2
                            p1(i,j) = pin(t);
                            %p2(i,j) = pin(t);
                        elseif j >= jd_c && j <= jd_c2
                            p1(i,j) = pin(t);
                            %p2(i,j) = pin(t);
                        elseif j >= jd_d && j <= jd_d2
                            p1(i,j) = pin(t);
                            %p2(i,j) = pin(t);
                        elseif j >= jd_e && j <= jd_e2
                            p1(i,j) = pin(t);
                            %p2(i,j) = pin(t);
                        else
                            p1(i,j) = cp1 .* p2(i,j) - cp2 .* (u2(i,j) - u2(i-1,j) + v2(i,j) - v2 (i,j-1));
                        end
                    else
                        p1(i,j) = cp1 .* p2(i,j) - cp2 .* (u2(i,j) - u2(i-1,j) + v2(i,j) - v2 (i,j-1));
                    end 
                else
                    p1(i,j) = cp1 .* p2(i,j) - cp2 .* (u2(i,j) - u2(i-1,j) + v2(i,j) - v2 (i,j-1));
                end
            end
        end
    end
    i = 1;
    for j = 2 : jx
        p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
%         p1(i + 1, j) = p1_taihi(i + 1, j);
    end
    i = ix + 1;
    for j = 2 : jx
        p1(i,j) = a .* (p1(i - 1,j) - p2(i,j)) - b .* (p1(i-2,j) - 2*p2(i-1,j) + p3(i,j)) - c .* (p2(i - 2,j) - p3(i-1,j)) + d .* p2(i -1,j) - e .* p3(i-2,j);
%          p1(i - 1, j) = p1_taihi(i - 1, j);
    end
    j = 1;
    for i = 2 : ix
        p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
%         p1(i ,j + 1) = p1_taihi(i, j + 1);
    end
    j = jx + 1;
    for i = 2 : ix
        p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
%         p1(i, j - 1) = p1_taihi(i,j - 1);
    end
    
    switch(boundary)
        case 0
            i = 1;
                for j = 2 : jx
                    p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
%                   p1(i + 1, j) = p1_taihi(i + 1, j);
                end
            i = ix + 1;
                for j = 2 : jx
                    p1(i,j) = a .* (p1(i - 1,j) - p2(i,j)) - b .* (p1(i-2,j) - 2*p2(i-1,j) + p3(i,j)) - c .* (p2(i - 2,j) - p3(i-1,j)) + d .* p2(i -1,j) - e .* p3(i-2,j);
                    p1(i - 1, j) = p3(i-1,j);
                end
            j = 1;
                for i = 2 : ix
                    p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
%                   p1(i ,j + 1) = p1_taihi(i, j + 1);
                end
            j = jx + 1;
                for i = 2 : ix
                    p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
%                   p1(i, j - 1) = p1_taihi(i,j - 1);
                end
        case 1
            i = 1;
                for j = 2 : jx
                    p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
                    p1(i + 1, j) = p3(i + 1, j);
                end
            i = ix + 1;
                for j = 2 : jx
                    p1(i,j) = a .* (p1(i - 1,j) - p2(i,j)) - b .* (p1(i-2,j) - 2*p2(i-1,j) + p3(i,j)) - c .* (p2(i - 2,j) - p3(i-1,j)) + d .* p2(i -1,j) - e .* p3(i-2,j);
                    p1(i - 1, j) = p3(i - 1, j);
                end
            j = 1;
                for i = 2 : ix
                    p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
                    p1(i ,j + 1) = p3(i, j + 1);
                end
            j = jx + 1;
                for i = 2 : ix
                    p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
                    p1(i, j - 1) = p3(i,j - 1);
                end
        case 2
            i = 1;
            for j = 2 : jx
                    p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
%                   p1(i + 1, j) = p1_taihi(i + 1, j);
                end
            i = ix + 1;
                for j = 2 : jx
                    p1(i,j) = a .* (p1(i - 1,j) - p2(i,j)) - b .* (p1(i-2,j) - 2*p2(i-1,j) + p3(i,j)) - c .* (p2(i - 2,j) - p3(i-1,j)) + d .* p2(i -1,j) - e .* p3(i-2,j);
%                    p1(i - 1, j) = p1_taihi(i - 1, j);
                end
            j = 1;
                for i = 2 : ix
                    p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
%                   p1(i ,j + 1) = p1_taihi(i, j + 1);
                end
            j = jx + 1;
                for i = 2 : ix
                    p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
%                   p1(i, j - 1) = p1_taihi(i,j - 1);
                end
    end
    if hekomi == 1
        for i = b_x : b_x + h_x
            for j = 1 : w_x
                p1(i,j) = 0;
                u1(i,j) = 0;
                v1(i,j) = 0;
            end
        end
        i = b_x - 1;
        for j = 1 : w_x + 1
            p1(i,j) = p3(i,j);
        end
        j = w_x + 1;
        for i = b_x - 1 : b_x + h_x + 1
            p1(i,j) = p3(i,j);
        end
        i = b_x + h_x + 1;
        for j = 1 : w_x + 1
            p1(i,j) = p3(i,j);
        end
    end
    
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            p3(i,j) = p2(i,j);
            p2(i,j) = p1(i,j);
        end
    end
    for i = 1 : ix
        for j = 2 : jx
            u1(i,j) = cv1 .* u2(i,j) - cv2 .* (p2(i + 1, j) - p2(i,j));
        end
    end
    for i = 2 : ix
        for j = 1 : jx
            v1(i,j) = cv1 * v2(i,j) - cv2 * (p2(i, j + 1) - p2(i,j));
        end
    end
    
    switch(boundary)
        case 0
            i = 1;
                for j = 2 : jx
%                     u1(i,j) = 0;
%                     v1(i,j) = 0;
                end
            i = ix + 1;
                for j = 2 : jx
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            j = 1;
                for i = 2 : ix
%                     u1(i,j) = 0;
%                     v1(i,j) = 0;
                end
            j = jx + 1;
                for i = 2 : ix
%                     u1(i,j) = 0;
%                     v1(i,j) = 0;
                end
        case 1
            i = 1;
                for j = 2 : jx
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            i = ix + 1;
                for j = 2 : jx
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            j = 1;
                for i = 2 : ix
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            j = jx + 1;
                for i = 2 : ix
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
    end
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            u2(i,j) = u1(i,j);
            v2(i,j) = v1(i,j);
        end
    end
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            pressure(i,j) = (p1(i,j) .^2).^0.5 ;
        end
    end
    
    p_keisoku_taihi(t) = p1(10,10);
    disp(t);
    time
    y = 1 : jx + 1 ;
    x = 1 : ix + 1;
    x_p = x * dx;
    y_p = y * dx;
    if mode_plot == 0
        if mod(t,10) == 0
        imagesc(y_p,x_p,pressure(x,y));
            colorbar
            %colormap gray ;
%             title(['pressure when ',num2str(time),'seconds have passed'])
            title(['pressure when frequency =',num2str(freq),'Hz&&', num2str(time), '(s) have passed'])
            xlabel('y(mm)')
            ylabel('x(mm)')
        grid on;
        drawnow
        end
    end
    if mode_plot == 1
        if t == tx
            y_half = round(jx / 2);
            plot(x_p,p1(x, y_half));
             min_v = min(p1(x, y_half));
             max_v = max(p1(x, y_half));
             for i = 1 : ix + 1
                 if p1(i,y_half) == min_v
                     min_x = i * dx;
                 elseif p1(i,y_half) == max_v
                     max_x = i * dx;
                 end
             end
            hacho = abs(max_x - min_x) * 2;
            min_x
            max_x
            hacho
            min_v
            v_c = hacho * freq;
            v_c
            crn
            break;
        end
    end
    if mode_plot ==2
        if mod(t,10) == 0
        plot(x_p,pressure(x, 10));
        grid on;
        drawnow
        end
    end
    if mode_plot == 3
        if t == 5000
            t_x = 1 : t;
            time = t_x * dt;
%            plot(time,p_keisoku_taihi(t_x));
%             bekutoru(time) = p_keisoku_taihi(t_x);
%             X = bekutoru(time);
            plot(time,p_keisoku_taihi(t_x));
            Fs = 1000;            % Sampling frequency                    
            T = 1/Fs;             % Sampling period       
            L = 1500;             % Length of signal
            ft = (0:L-1)*T;        % Time vector
            Y = fft(X);
%             plot(1:5000,Y);
            %P2 = abs(Y/L);
            %P1 = P2(1:L/2+1);
            %P1(2:end-1) = 2*P1(2:end-1);
            %f = Fs*(0:(L/2))/L;
            %plot(f,P1);
            %title('Single-Sided Amplitude Spectrum of X(t)')
            %xlabel('f (Hz)')
            %ylabel('|P1(f)|')
            break;
        end
    end
end
end
if mode == 1
    if range ~= 2
        for i = id : id2
            for j = jd : jd_2
                p1(i,j) = 100;
            end
        end
    else
        for i = id : id2
            for j = jd_a : jd_a2
                p1(i,j) = 100;
            end
            for j = jd_b : jd_b2
                p1(i,j) = 100;
            end
            for j = jd_c : jd_c2
                p1(i,j) = 100;
            end
            for j = jd_d : jd_d2
                p1(i,j) = 100;
            end
            for j = jd_e : jd_e2
                p1(i,j) = 100;
            end
        end
    end
    if range == 0
        for i = b_x : b_x + h_x
            for j = 1 : w_x
                p1(i,j) = 50;
            end
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
if mode == 2
    disp(w); 
    disp(freq);
    x = 0 : 0.1 : 2 * 3.1415;
    y = sin(x);
    Y = fft(y);
    plot(Y)
    crn
end
if mode == 3
  for freq = 100 : 100 : 1000000
    theta = pi:0.01:3*pi;
    %dhi = (sin(5*3.14*b*freq*sin(theta)/340))/(5 * sin(3.14*b*freq*sin(theta))/340)
    dhi1 = (sin(5*3.14*b_haba*freq*sin(theta)/340)) ./(5 * sin(3.14*b_haba*freq*sin(theta)/340)) ;
    if mod(freq,5000) == 0
        polarplot(theta,abs(dhi1));
        title(['Directivity when frequency is',num2str(freq),'Hz'])
        drawnow
        pause(1);
    end
  end
end
if mode == 4
    hasu = 2 * pai * freq / 340;
    x = 1 : ix + 1;
    x_p = x * dx;
    t_time = 0.01;
    j = sqrt(-1);
    p = exp(j *(hasu * x_p - w * t_time));
    plot(x_p,p)
end