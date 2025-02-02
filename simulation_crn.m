range = 0;
boundary = 0;
mode = 0; %0...通常シミュレーション 1...初期印加位置表示 2...変数表示3...指向性の極グラフ
mode_plot = 6; %プロットモード選択 0...カラーマップ進行 1...xプロット 2...xプロット進行 3...ある地点の時間変化 4..先行研究 5...ある地点のパワースペクトル
% 6...3を細かい時間で追う
hekomi = 0;
sweep = 1;
SC = 0;%励振関数 ０なら連続1ならガウシアン2ハニング3正弦波数波4スイープ

f1 = figure;
% f2 = figure;

c = 340; %音速
c0 = 340;
% rou0 = 1.293; %密度（kg/m^3
rou0 = 1.293;
% rou0 = 1000;?
freq_param = 1;
freq_a = 1000;
freq_start = 000*freq_param;
freq_add = 4000*freq_param;
freq = freq_a*freq_param;
freq_abs = freq_a*freq_param;
ramuda = c0 / freq;
dx_param = 0.005; %0.05-0.025 
dx = ramuda*dx_param; % λの20-30分の一
% dt = dx / (5 * c);
% dt = dx*0.15 / (c0);
crn_param = 0.2;
dt = dx*crn_param/ (c0);
% クー数から条件を立てる]

cal_time = 0.18;
if range == 0
    xrange = 2.01;
%    xrange = 0.52;
    yrange = 0.02;
%     yd = 0.01 - 2*dx;
%     yd2 = 0.01 + 2*dx;
    yd = dx;
    yd2 = 0.02 - dx;
    xd = 3 * dx;
    xd2 = 4 * dx;
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
absp0 = - 0.5; % 吸収係数
b_po = 0.5 ; %凹み位置
h = 0.05;%凹み幅
w = 0.015;%凹みふかさ


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
WN = 1;
W_end = round(2*WN/(freq*dt)) - 1 ;
pin = zeros(1,tx+100);
p_keisoku_taihi = zeros(1,tx);
crn =(c0 * dt)/dx ; %クーラン数
dd = 199/200 ;%higdons absorption boundary
pai = 3.1415 ;

% 編集点↓
gensui0 = (freq*absp0) / (8.686*c0); % 減衰係数
hasu_o0 = 2*pai*freq/c0 ;%実際の波数
hasu0 = sqrt(hasu_o0*hasu_o0 - gensui0^2);%損失ある場合の波数
c_m0 = 2*pai*freq/hasu0;%損失のある場合の音速
alpha0 = 2*hasu_o0*rou0*c_m0/hasu0;%吸収項
kap0 = c_m0^2*rou0;
cp1 = 1;
cp2 = rou0*c_m0^2*dt/dx;
%cp2 = crn;
cv1 = 1;
cv2 = dt / (rou0 * dx);
% 編集点↑

a = 2*(crn - 1)/(crn + 1);
b = ((crn - 1)/(crn + 1))*2;%要検討要検討
c = 2 * (crn - 1) / (crn + 1) * dd;
d = dd * 2;
e = dd *2; %要検討要検討
a1 = 1 / cos(0);
a2 = 1 / cos(0);
d1 = 0.005;
d2 = 0.005;

% 編集点↓
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
% 編集点↑

y_half = round(jx / 2);
p_kei = zeros(1, tx);
% y_half_g = round(y_half / dx);
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
w = 2 * pai * freq ;
    switch SC
        case 0
            for tc = 1 : tx+100 
                time = tc * dt;
                tr = w * dt * tc;
                    if sweep == 1
                        %freq = 2000 + 7000*(time/cal_time);
                        freq = freq_start + freq_add*(time/cal_time);
                    end
%                 freq = 5000;
                w = 2 * pai * freq ;
                j = sqrt(-1);
                if tc < pai / (w * dt)
                     pin(tc) = sin(w * time);
%                     pin(tc) = exp(j * w * time);
                else
                     pin(tc) = sin(w * time);
%                     pin(tc) = exp(j * w * time);
%                   pin(t) = 0;
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
        case 4
            for tc = 1 : tx+100 
                time = tc * dt;
                tr = w * dt * tc;
                    if sweep == 1
                        %freq = 2000 + 7000*(time/cal_time);
                        freq = freq_start + freq_add*(time/cal_time);
                    else
                        freq = freq_abs;
                    end
%                 freq = 5000;
                w = 2 * pai * freq ;
                j = sqrt(-1);
                if tc < pai / (w * dt)
%                     pin(tc) = sin(w * time);
                    pin(tc) = exp(j * w * time - j * pi/2);
                else
%                     pin(tc) = sin(w * time);
                    pin(tc) = exp(j * w * time - j * pi / 2);
%                   pin(t) = 0;
                end
            end
    end

if mode == 0
for t = 1: tx
   time = t * dt;
   if sweep == 1
       freq = 2000 + 7000*(time/cal_time);
   end
   if crn >= 1
       disp("クーラン数が不適切です。");
       break;
   end
   
	ramuda_2 = c0 / freq;
    dx_2 = ramuda/25;
    dt_2 = dx_2 / (5 * c0);
    gensui0 = (freq*absp0) / (8.686*c0); % 減衰係数
    hasu_o0 = 2*pai*freq/c0 ;%実際の波数
    hasu0 = sqrt(hasu_o0*hasu_o0 - gensui0^2);%損失ある場合の波数
    c_m0 = 2*pai*freq/hasu0;%損失のある場合の音速
    alpha0 = 2*hasu_o0*rou0*c_m0/hasu0;%吸収項
    kap0 = c_m0^2*rou0;
    cp1 = 1;
    cp2 = rou0*c_m0^2*dt_2/dx_2;
%     cv1 = 1;
%     cv2 = dt / (rou0 * dx);
    cv1 = (2*rou0-alpha0*dt_2)/(2*rou0+alpha0*dt_2);
    cv2 = (2*dt_2)/((2*rou0+alpha0*dt_2)*dx_2);
    ca0 = (c_m0 * dt_2 - dx_2) / ((c_m0 * dt_2 + dx_2));
    ca1 = (dx_2 * 2) / (c_m0 * dt_2 + dx_2);
    ca2 = (dx_2 * c_m0 * dt_2 * c_m0 * dt_2) / (dx_2 * dx_2 * 2 * (c_m0 * dt_2 + dx_2));
    cah1 = (a1 * c_m0 * dt_2 - dx_2) / (a1 * c_m0 * dt_2 + dx_2);
    cah2 = (a2 * c_m0 * dt_2 - dx_2) / (a2 * c_m0 * dt_2 + dx_2);
    a = cah1 + cah2;
    b = cah1 * cah2;
    c = cah1 * (1 - d2) + cah2 * (1 - d1);
    d = ((1 - d1) + (1 - d2));
    e = ((1 - d1) * (1 - d2));

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
                        s = 7;
                        if j >= jd_a && j <= jd_a2
                            p1(i,j) = pin(t);
                            %p2(i,j) = pin(t);
                        elseif j >= jd_b && j <= jd_b2
                            p1(i,j) = pin(t + s);
                            %p2(i,j) = pin(t);
                        elseif j >= jd_c && j <= jd_c2
                            p1(i,j) = pin(t + 2*s);
                            %p2(i,j) = pin(t);
                        elseif j >= jd_d && j <= jd_d2
                            p1(i,j) = pin(t + 3*s);
                            %p2(i,j) = pin(t);
                        elseif j >= jd_e && j <= jd_e2
                            p1(i,j) = pin(t + 4*s);
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
                    %p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
                    p1(i + 1, j) = p2(i + 1, j);
                    p1(i,j) = 0;
                end
            i = ix + 1;
                for j = 2 : jx
                    %p1(i,j) = a .* (p1(i - 1,j) - p2(i,j)) - b .* (p1(i-2,j) - 2*p2(i-1,j) + p3(i,j)) - c .* (p2(i - 2,j) - p3(i-1,j)) + d .* p2(i -1,j) - e .* p3(i-2,j);
                    p1(i - 1, j) = p2(i - 1, j);
                    p1(i,j) = 0;
                end
            j = 1;
                for i = 2 : ix
                    %p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
                    p1(i ,j + 1) = p2(i, j + 1);
                    p1(i,j) = 0;
                end
            j = jx + 1;
                for i = 2 : ix
                    %p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
                    p1(i, j - 1) = p2(i,j - 1);
                    p1(i,j) = 0;
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
%                 p1(i,j) = 0;
                u1(i,j) = 0;
                v1(i,j) = 0;
            end
        end
        i = b_x - 1;
        for j = 1 : w_x + 1
%             p1(i,j) = p2(i,j);
        end
        j = w_x + 1;
        for i = b_x - 1 : b_x + h_x + 1
%             p1(i,j) = p2(i,j);
        end
        i = b_x + h_x + 1;
        for j = 1 : w_x + 1
%             p1(i,j) = p2(i,j);
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
            pressure(i,j) = (real(p1(i,j)) .^2).^0.5 ;
        end
    end
    if mod(t,5)==0
%         計測点
        p_keisoku_taihi(t/5) = p1(5, y_half);
    end
    if mod(t,500) == 0
        disp(t);
        time
    end
    
    y = 1 : jx + 1 ;
    x = 1 : ix + 1;
    x_p = x * dx;
    y_p = y * dx;
    if mode_plot == 0
        if mod(t,100) == 0
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
%         disp(p1(4,y_half));
        if mod(t,50) == 0
            plot(x_p(1:0.5/dx),p1(x(1:0.5/dx),y_half))
            grid on;
            drawnow
        end
        %if t > 3500
        if time > 4.5e-04
            if real(p1(1,y_half)) < 0.001 && real(p1(4,y_half)) > -0.001
            %plot(x,p1(x_p, y_half));
            plot(x_p(1:0.5/dx),p1(x(1:0.5/dx), y_half));
%              min_v = min(p1(x, y_half))
            min_v = 0;
            for i = 1 : ix + 1
                if real(p1(i,y_half)) < min_v
                    min_v = real(p1(i,y_half));
                end
            end
             max_v = max(real(p1(x, y_half)))
             min_v
             for i = 1 : ix + 1
                 difference_min(i) = abs(real(p1(i,y_half)) - min_v);
                 difference_max(i) = abs(real(p1(i,y_half)) - max_v);
             end
             for i = 1 : ix + 1
                 if difference_min(i) == min(difference_min)
                     minor_i = i;
                 end
                 if difference_max(i) == min(difference_max)
                     major_i = i;
                 end
             end
             min_x = minor_i * dx;
             max_x = major_i * dx;
%              for i = 1 : ix + 1
%                  if real(p1(i,y_half)) == min_v
%                      min_x = i * dx;
%                  elseif real(p1(i,y_half)) == max_v
%                      max_x = i * dx;
%                  end
%              end
            
                hacho = abs(max_x - min_x) * 2;
                min_x
                max_x
                hacho
                min_v
                v_c = hacho * freq_a
                break;
            end
        end
    end
    if mode_plot ==2
        if mod(t,10) == 0
%         plot(x_p,pressure(x, y_half));
%         plot(x_p,p1(x, y_half));
        plot(x_p(1:0.5/dx),p1(x(1:0.5/dx),y_half))
        grid on;
        drawnow
        end
    end
    if mode_plot == 3
        if t == 5000
%         if mod(t,50) == 0
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
    if mode_plot == 4
        p_kei(t) = p1(1, y_half);
        p_kei(t)
        if time == cal_time - 0.001
            t_x = 1 : tx;
            plot(t_x * dt , p_kei(t_x))
            disp("計算が終了しました")
        end
    end
    if mode_plot == 5
        if t == tx
            t_x = 0 : 10: t;
            time = t_x * dt;
            p_keisoku_spec = abs(p_keisoku_taihi).^2;
            plot(time,p_keisoku_spec(t_x));
            Fs = 1000;            % Sampling frequency                    
            T = 1/Fs;             % Sampling period       
            L = 1500;             % Length of signal
            ft = (0:L-1)*T;        % Time vector
            Y = fft(p_keisoku_spec,L);
            disp(size(Y))
            plot(1:L,Y)
            break;
        end
    end
    if mode_plot == 6
        if mod(t,500) == 0
            figure(f1);
            t_x = 1 : 1: t/10;
            time = t_x * dt * 10;
           p_keisoku_spec = (abs(p_keisoku_taihi)).^2;
           plot(time, p_keisoku_taihi(t_x)); 
           grid on;
           drawnow;
        end
%        if time > 0.01
%            if mod(t,10) == 0
%                p_keisoku_spec(t/10)
%                time
%                break;
%            end
%        end
%             if t == tx - rem(tx,10)
            if t == tx - rem(tx,10)
                disp(size(time));
                disp(size(p_keisoku_spec));
                disp("終了")
                %csv_array = [time; p_keisoku_spec];
                p_keisoku_spec_col = p_keisoku_spec.';
                p_keisoku_taihi = p_keisoku_taihi.';
                %dlmwrite('kairyouheko500sweep2to7.csv', p_keisoku_spec_col, 'precision', '%.10f', 'delimiter', ',')
                dlmwrite('pownoload_200.csv', p_keisoku_taihi, 'precision', '%.10f', 'delimiter', ',')
                break;
            end
    end
    if mode_plot == 7
        if time > 0.001
            plot(x_p,p1(x, y_half));
            break;
        end
    end
    if mode_plot == 8
        sokutei_point = 300;
        if abs(p_keisoku_taihi(sokutei_point)) > 0
            sokutei_point * dx
            disp(sokutei_point * dx / time);
            break;
        end
        if mod(t,5) == 0
%         plot(x_p,pressure(x, y_half));
            plot(x_p,p1(x, y_half));
            grid on;
            drawnow
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
    ramuda
    dx
    dt
end
if mode == 3
  for freq = 100 : 100 : 10000000
    theta = pi:0.01:3*pi;
    dhi = 0.05;
    enu = 4;
    b_haba = dhi / enu;
    %dhi = (sin(5*3.14*b*freq*sin(theta)/340))/(5 * sin(3.14*b*freq*sin(theta))/340)
    dhi1 = (sin(enu*3.14*b_haba*freq*sin(theta)/340)) ./(enu * sin(3.14*b_haba*freq*sin(theta)/340)) ;
    if mod(freq,1000) == 0
        polarplot(theta,abs(dhi1));
        title(['Directivity when frequency is',num2str(freq),'Hz'])
        drawnow
        pause(0.7);
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