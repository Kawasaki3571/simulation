range = 0;
boundary = 0;
mode = 0; %0...通常シミュレーション 1...初期印加位置表示 2...変数表示3...指向性の極グラフ
mode_plot = 0; %プロットモード選択 0...カラーマップ進行 1...xプロット 2...xプロット進行 3...ある地点の時間変化 4..先行研究 5...ある地点のパワースペクトル
% 6...3を細かい時間で追う
hekomi =  0;

sweep = 0;
SC = 0;%励振関数 ０なら連続1ならガウシアン2ハニング3正弦波数波4スイープ5インパルス
stair = 3;
circlecase = 1;

f1 = figure;
% f2 = figure;

c = 340; %音速
c0 = 340;
% rou0 = 1.293; %密度（kg/m^3
rou0 = 1.293;
% rou0 = 1000;
freq_param = 1;
freq_a = 10000;
freq_start = 000*freq_param;
freq_add = 12000*freq_param;
% freq = freq_a*freq_param;

freq = freq_a;
freq = 30000;

freq_abs = freq_a*freq_param;
ramuda = c0 / freq;
dx_param = 0.005; %0.05-0.025

% dx_param = 0.01;

dx = ramuda*dx_param; % λの20-30分の一
dx = 0.001 ;

crn_param = 0.2;
crn = crn_param;
dt = dx*crn_param/ (c0);
% クー数から条件を立てる]

cal_time = 0.18;


if range == 0
    xrange = 0.41;
%    xrange = 0.51;
    yrange = 0.05;
    zrange = 0.05;
    yd = yrange / 2;
    yd2 = yd;
    zd = zrange / 2;
    zd2 = zd;
    xd = 3 * dx;
    xd2 = 4 * dx;
end
    
absp0 = - 0.5; % 吸収係数
b_po = 0.3 ; %凹み位置
h = 0.005;%凹み幅
w = 0.005;%凹みふかさ
b_de = 0.02;

ix = round(xrange / dx); %x空間感覚の数
jx = round(yrange / dx); %y空間感覚の数
kx = round(zrange / dx);
tx = fix(cal_time / dt ); %時間感覚の数
b_x = round(b_po / dx) + 1;
h_x = round(h / dx) + 1;
w_x = round(w / dx) + 1;
bd_x = round(b_de / dx) + 1;
%  １おわ
td = 15; % 周波数の代入（？）
id = round(xd / dx);%xの位置（？）
id2 = round(xd2 / dx ) ;
%sokuteiten_x_g = round(sokuteiten_x / dx);
%sokuteiten_y_g = round(sokuteiten_y / dx);
jd = round(yd / dx) + 1  ; %y空間感覚の代入（？）
jd_2 = round(yd2/dx) + 1;
kd = round(zd / dx);
kd_2 = round(zd2 / dx);

%ix1 = round(x1/dx);
%ix2 = round(x2/dx);
%jx1 = round(y1/dx);
%jx2 = round(y2/dx);

p1 = zeros(ix+1,jx+1,kx+1);
p2 = zeros(ix+1,jx+1,kx+1);
p3 = zeros(ix+1,jx+1,kx+1);
u1 = zeros(ix+1,jx+1,kx+1);
u2 = zeros(ix+1,jx+1,kx+1);
v1 = zeros(ix+1,jx+1,kx+1);
v2 = zeros(ix+1,jx+1,kx+1);
w1 = zeros(ix+1,jx+1,kx+1);
w2 = zeros(ix+1,jx+1,kx+1);
pressure = zeros(ix + 1 , jx + 1,kx + 1);
WN = 1;
W_end = round(2*WN/(freq*dt)) - 1 ;
pin = zeros(1,tx+100);
p_keisoku_taihi = zeros(1,round(tx / 5));

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
z_half = round(kx / 2);
p_kei = zeros(1, tx);
% y_half_g = round(y_half / dx);
for i = 1 : ix + 1
    for j = 1 : jx + 1
        for k = 1 : kx + 1;
            p1(i,j,k) = 0;
            p2(i,j,k) = 0;
            p3(i,j,k) = 0;
            u1(i,j,k) = 0;
            u2(i,j,k) = 0;
            v1(i,j,k) = 0;
            v2(i,j,k) = 0;
            w1(i,j,k) = 0;
            w2(i,j,k) = 0;
            p1_taihi(i,j,k) = 0;
        end
    end
end
x=0:dx:xrange; 
y=0:dx:yrange;
z=0:dx:zrange;
%%timeloop
if xd > xd2
    error("xrangeがおかしい")
end
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
                %freq = 10000;
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
                    pin(tc) = exp(j * w * time);
                else
%                     pin(tc) = sin(w * time);
                    pin(tc) = exp(j * w * time);
%                   pin(t) = 0;
                end
            end
        case 5
            pulse_second = 0.0001;
            for tc = 1 : tx + 100
                time = tc * dt;
                if time < pulse_second
                    pin(tc) = 1;
                else
                    pin(tc) = 0;
                end
            end
    end

if mode == 0
for t = 1: tx
   time = t * dt;
   if sweep == 1
%       freq = 2000 + 7000*(time/cal_time);
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

    for i = 2 : ix + 1
        for j = 2 : jx + 1
            for k = 2 : kx + 1
                if range ~= 2
                    if t <= tx && i >= id && i <= id2 && j >= jd && j <= jd_2 && k >= kd && k <= kd_2
                        p1(i,j,k) = pin(t);
                        p2(i,j,k) = pin(t);
                    else
                        p1(i,j,k) = cp1 .* p2(i,j,k) - cp2 .* (u2(i,j,k) - u2(i-1,j,k) + v2(i,j,k) - v2 (i,j-1,k) + w2(i,j,k) - w2 (i,j,k-1));
                    end
                else
                    if t <= tx
                        if i >= id && i <= id2
                        else
                            p1(i,j) = cp1 .* p2(i,j) - cp2 .* (u2(i,j) - u2(i-1,j) + v2(i,j) - v2 (i,j-1));
                        end 
                    else
                        p1(i,j) = cp1 .* p2(i,j) - cp2 .* (u2(i,j) - u2(i-1,j) + v2(i,j) - v2 (i,j-1));
                    end
                end
            end
        end
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
                    u1(i, j) = 0;
                    v1(i, j) = 0;
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
            for i = 1 : 1
                for j = 1 : jx + 1
                    %p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            end
            i = ix + 1;
                for j = 1 : jx + 1
                    %p1(i,j) = a .* (p1(i - 1,j) - p2(i,j)) - b .* (p1(i-2,j) - 2*p2(i-1,j) + p3(i,j)) - c .* (p2(i - 2,j) - p3(i-1,j)) + d .* p2(i -1,j) - e .* p3(i-2,j);
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            for j = 1 : 1
                for i = 1 : ix + 1
                    %p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            end
            j = jx + 1 ;
                for i = 1 : ix + 1
                    %p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
        case 2
            i = 1;
                for j = 2 : jx 
                    p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
%                   p1(i + 1, j) = p1_taihi(i + 1, j);
                end
            i = ix + 1;
                for j = 1 : jx
                    p1(i,j) = a .* (p1(i - 1,j) - p2(i,j)) - b .* (p1(i-2,j) - 2*p2(i-1,j) + p3(i,j)) - c .* (p2(i - 2,j) - p3(i-1,j)) + d .* p2(i -1,j) - e .* p3(i-2,j);
%                    p1(i - 1, j) = p1_taihi(i - 1, j);
                end
            j = 1;
                for i = 1 : ix
                    p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
%                   p1(i ,j + 1) = p1_taihi(i, j + 1);
                end
            j = jx + 1;
                for i = 2 : ix 
                    p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
%                   p1(i, j - 1) = p1_taihi(i,j - 1);
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
        for j = 1 : jx + 1
            for k = 1 : kx + 1
                u1(i,j,k) = cv1 .* u2(i,j,k) - cv2 .* (p2(i + 1, j, k) - p2(i,j,k));
            end
        end
    end
    for i = 1 : ix + 1
        for j = 1 : jx
            for k = 1 : kx + 1
                v1(i,j,k) = cv1 .* v2(i,j,k) - cv2 .* (p2(i, j+1, k) - p2(i,j,k));
            end
        end
    end
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            for k = 1 : kx
                w1(i,j,k) = cv1 .* w2(i,j,k) - cv2 .* (p2(i, j, k+1) - p2(i,j,k));
            end
        end
    end
    if hekomi == 1
        for i = b_x : b_x + h_x
            for j = 1 : w_x
                for k = 1 : bd_x
                    u1(i,j,k) = 0;
                    v1(i,j,k) = 0;
                    w1(i,j,k) = 0;
                end
            end
        end
    end
        
    if stair == 2
        center_x = 0.1;
        center_y = 0.1;
        kei = 0.1;
        for x1_s = 0 : dx : xrange
            y1_s = center_y + sqrt(kei^2 - (x1_s - center_x)^2);
            i1 = round(x1_s/dx) + 1;
            j1 = round(y1_s/dx) + 1;
            u1(i1, j1) = 0;
            v1(i1, j1) = 0;
            for jstair = j1 : jx + 1
                u1(i1, jstair) = 0;
                v1(i1, jstair) = 0;
            end
        end
        for x1_s = 0 : dx : xrange
            y1_s = center_y - sqrt(kei^2 - (x1_s - center_x)^2);
            i1 = round(x1_s/dx) + 1;
            j1 = round(y1_s/dx) + 1;
            u1(i1, j1) = 0;
            v1(i1, j1) = 0;
            for jstair = 1 : j1
                u1(i1, jstair) = 0;
                v1(i1, jstair) = 0;
            end
        end
    end
    if stair == 3
        for istair = 1 : 2
            for jstair = 1 : jx + 1
                for kstair = 1 : kx + 1
                    u1(istair, jstair, kstair) = 0;
                    v1(istair, jstair, kstair) = 0;
                end
            end
        end
        for istair = 3 : ix - 1
            for jstair = 1 : 2
                for kstair = 1 : kx + 1
                    u1(istair, jstair, kstair) = 0;
                    v1(istair, jstair, kstair) = 0;
                end
            end
            for jstair = 3 : jx - 1
                for kstair = [1 2 kx kx+1]
                    u1(istair, jstair, kstair) = 0;
                    v1(istair, jstair, kstair) = 0;
                end
            end
            for jstair = jx  : jx + 1
                for kstair = 1 : kx + 1
                    u1(istair, jstair, kstair) = 0;
                    v1(istair, jstair, kstair) = 0;
                end
            end
        end
        for istair = ix : ix + 1
            for jstair = 1 : jx + 1
                for kstair = 1 : kx + 1
                    u1(istair, jstair, kstair) = 0;
                    v1(istair, jstair, kstair) = 0;
                end
            end
        end
        if circlecase == 1
            for istair = 1 : ix + 1
                center_z = zrange/2;
                center_y = yrange/2;
                kei = zrange / 2;
                for z1_s = 0 : dx : zrange
                    y1_s = center_y + sqrt(kei^2 - (z1_s - center_z)^2);
                    i1 = round(z1_s/dx) + 1;
                    j1 = round(y1_s/dx) + 1;
                    %u1(istair, i1, j1) = 0;
                    %v1(istair, i1, j1) = 0;
                        for jstair = j1 : jx + 1
                            u1(istair, i1, jstair) = 0;
                            v1(istair, i1, jstair) = 0;
                            %v1(i1, jstair) = 0;
                        end
                end
                for z1_s = 0 : dx : zrange
                    y1_s = center_y - sqrt(kei^2 - (z1_s - center_z)^2);
                    i1 = round(z1_s/dx) + 1;
                    j1 = round(y1_s/dx) + 1;
                    %p1(i1, j1) = 0;
                    %v1(i1, j1) = 0;
                    for jstair = 1 : j1
                        u1(istair, i1, jstair) = 0;
                        v1(istair, i1, jstair) = 0;
                    end
                end
            end
        end
    end
    
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            for k = 1 : kx + 1
                u2(i,j,k) = u1(i,j,k);
                v2(i,j,k) = v1(i,j,k);
            end
        end
    end
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            for k = 1 : kx + 1
                pressure(i,j,k) = (real(p1(i,j,k)) .^2).^0.5 ;
            end
        end
    end
    if mod(t,5)==0
%         計測点
        p_keisoku_taihi(t/5) = p1(6, y_half, z_half);
%        p_keisoku_taihi(t/5) = p1(sokuteiten_x_g, sokuteiten_y_g);
    end
    if mod(t,50) == 0
        disp(t);
        time
    end
    
    y = 1 : jx + 1 ;
    x = 1 : ix + 1;
    z = 1 : kx + 1;
    x_p = x * dx;
    y_p = y * dx;
    z_p = z * dx;
    if mode_plot == 0
        if mod(t,10) == 0
%         imagesc(y_p,x_p,pressure(x,y));
            z = z_half;
            
            %plot_image = pressure(z)';
            for x_plot = 1 : ix + 1
                for y_plot = 1 : jx + 1
                    plot_image(x_plot, y_plot) = pressure(x_plot, y_plot, z_half);
                end
            end
            
            imagesc(y_p, x_p, plot_image(x, y));
            
            %surf(x_p, y_p, abs(p1(x, y)));
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
    if mode_plot == 6
        if mod(t,500) == 0
            figure(f1);
            t_x = 1 : 1: t/5;
            time = t_x * dt * 5;
           p_keisoku_spec = (abs(p_keisoku_taihi)).^2;
           plot(time, p_keisoku_taihi(t_x));
           grid on;
           drawnow;
        end
            if t == tx - rem(tx,10)
                disp(size(time));
                disp(size(p_keisoku_spec));
                disp("終了")
                %csv_array = [time; p_keisoku_spec];
                p_keisoku_spec_col = p_keisoku_spec.';
                p_keisoku_taihi = p_keisoku_taihi.';
                dlmwrite('3dtsu_300_4cm_1.7mm_180ms.csv', p_keisoku_taihi, 'precision', '%.15f', 'delimiter', ',')
                %dlmwrite('pow500_0to4_1cm.csv', p_keisoku_taihi, 'precision', '%.10f', 'delimiter', ',')
                break;
            end
    end
    if mode_plot == 7
        if time > 0.001
            plot(x_p,p1(x, y_half));
            break;
        end
    end
    if mode_plot == 9
        if mod(t,50) == 0
            [X,Y,Z] = meshgrid(x,y,z);
            xslice = 1;   
            yslice = 1;
            zslice = 1;
            slice(X,Y,Z,p1(X,Y,Z),xslice,yslice,zslice)
            colorbar
            title(['3d'])
            xlabel('y(mm)')
            ylabel('x(mm)')
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
%        pause(0.7);
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
if mode == 5
    for x1 = 0 : dx : 0.1
        y1 = 0.1 - x1;
        i1 = round(x1/dx) + 1;
        j1 = round(y1/dx) + 1;
        p1(i1, j1) = 1;
    end
    for x1 = 0.1 : dx : 0.2
        y1 = x1 - 0.1;
        i1 = round(x1/dx) + 1;
        j1 = round(y1/dx) + 1;
        p1(i1, j1) = 0.5;
    end
    for x1 = 0.1 : dx : 0.2
        y1 = 0.3 - x1;
        i1 = round(x1/dx) + 1;
        j1 = round(y1/dx) + 1;
        p1(i1, j1) = -0.5;
    end
    for x1 = 0 : dx : 0.1
        y1 = x1 + 0.1;
        i1 = round(x1/dx) + 1;
        j1 = round(y1/dx) + 1;
        p1(i1, j1) = -1;
    end
    y = 1 : jx + 1 ;
    x = 1 : ix + 1;
    x_p = x * dx;
    y_p = y * dx;
    imagesc(y_p,x_p,p1(x,y));
    colorbar
    %colormap gray ;
    %   title(['pressure when ',num2str(time),'seconds have passed'])
    title(['pressure when frequency =',num2str(freq),'Hz&&', num2str(time), '(s) have passed'])
    xlabel('y(mm)')
    ylabel('x(mm)')
    grid on;
end
