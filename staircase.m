% 二次元における音波の伝播をシミュレーションするスクリプト

% 境界条件やプロット法を関数等にまとめておらず下の変数の値を操作し各条件を設定する。
% !!!他の使用者のための可読性向上のためリファクタリングが必要!!!

range = 1;%1...チューブ2...正方形5...ロボットハンド
boundary = 1; % 0...使わない1...四方完全反射境界2...四方吸収境界
mode = 0; %0...通常シミュレーション 1...初期印加位置表示 2...変数表示3...指向性の極グラフ6..へこみ初期印加位置表示
mode_plot = 0; %プロットモード選択 0...カラーマップ進行 1...xプロット 2...xプロット進行 3...ある地点の時間変化 4..先行研究 5...ある地点のパワースペクトル
% 6...3を細かい時間で追う、音圧をcsvに出力
hekomi = 0;%0でチューブモデルの凹みなし、1であり
sweep = 0;%0でスイープなし、1であり
SC = 0;%初期印加励振関数 ０なら連続1ならガウシアン2ハニング3正弦波数波4スイープ複素数5インパルス
stair = 1;%階段近似による多様な形状の境界実装0...ひし形1...円形3...四方完全反射境界5...ロボットハンド

f1 = figure;
% f2 = figure;

c = 340; %音速
c0 = 340;
rou0 = 1.293;
% rou0 = 1000;
freq_param = 1;
freq_a = 1000;
freq_start = 000*freq_param;%スイープの周波数
freq_add = 12000*freq_param;
% freq = freq_a*freq_param;

% 励振条件
ippa = 1;
% 特定時間における波形画像取得のための変数
snap1 = 0;
snap2 = 0;
snap3 = 0;
snap4 = 0;
snap_time1 = 0.075/c;
snap_time2 = 0.125/c;
snap_time3 = 0.175/c;
snap_time4 = 0.225/c;

freq = freq_a;

freq_abs = freq_a*freq_param;
ramuda = c0 / freq;
% dx_param = 0.002; %0.05-0.025
dx_param = 0.005;

dx = ramuda*dx_param; % λの20-30分の一が適正値とされている
dx = 0.001;

if range == 5%ロボットハンドなら領域が狭いためより細かく
    dx = 0.0005;
end
crn_param = 0.2;
dt = dx*crn_param/ (c0);

cal_time = 0.18;%シミュレーション時間
% cal_time = 0.36;
cal_time = 0.06;

%rangeごとのシミュレーション領域の定義
if range == 0
    xrange = 2.01;
    yrange = 0.02;
    yd = yrange / 2;
    yd2 = yd;
    xd = 3 * dx;
    xd2 = 4 * dx;
end
if range == 1
    xrange = 0.2;
    yrange = 0.2;
    xrange = 0.2*sqrt(2);
    yrange = 0.2*sqrt(2);
    xd = xrange / 2 ;
    xd2 = xrange / 2 ;
    yd = yrange / 2 ;
    yd2 = yrange / 2 ;
    sokuteiten_x = 0.1 * (3*sqrt(2)/4);
    sokuteiten_y = 0.1 * (3*sqrt(2)/4);
    sokuteiten_x = 0.15;
    sokuteiten_y = 0.1;
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
if range == 4
    xrange = .8;
    yrange = 0.02;
    xd = xrange / 2 ;
    xd2 = xrange / 2 ;
    
    xd = 100*dx;
    xd2 = xd;
    
    yd = yrange / 2 ;
    yd2 = yrange / 2 ;
    sokuteiten_x = 0.1 * (3*sqrt(2)/4);
    sokuteiten_y = 0.1 * (3*sqrt(2)/4);
    sokuteiten_x = 0.15;
    sokuteiten_y = 0.1;
end
if range == 5
    xrange = 0.06;
    yrange = 0.11;
    xd = 0.005 + 3*dx;
    xd2 = xd;
    yd = 3*dx;
    yd2 = yd;
end

x1 = 0.03;
x2 = 0.04;
y1 = 0.1;
y2 = 0.44;
t1 = 0;
t2 = 0;
speed = 0;
absp0 = - 0.5; % 吸収係数
b_po = 0.6 ; %凹み位置
h = 0.005;%凹み幅
w = 0.001;%凹みふかさ
b_po2 = 0.6 ; %凹み位置
h2 = 0.02;%凹み幅
w2 = 0.015;%凹みふかさ

ix = round(xrange / dx); %x空間感覚の数
jx = round(yrange / dx); %y空間感覚の数
tx = fix(cal_time / dt ); %時間感覚の数
b_x = round(b_po / dx) + 1;%
h_x = round(h / dx) + 1;
w_x = round(w / dx) + 2;
b_x_2 = round(b_po2 / dx) + 1;
h_x_2 = round(h2 / dx) + 1;
w_x_2 = round(w2 / dx) + 1;
id = round(xd / dx);%xの位置（？）
id2 = round(xd2 / dx ) ;
if range ~= 2
    jd = round(yd / dx)  ; %y空間感覚の代入（？）
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
p_keisoku_taihi = zeros(1,round(tx / 5));
p_keisoku_detail = zeros(1,round(tx/2));
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

% 音圧更新粒子速度更新吸収境界のための必要変数
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

y_half = round(jx / 2);
p_kei = zeros(1, tx);
% 音圧粒子速度行列初期化（matlabならなくても動くかも）
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
sn = 0;


% シミュレーションが動くための必要条件
if xd > xd2
    error("xrangeがおかしい")
end
if range ~= 2 && range ~= 5
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
% 初期印加波の定義
    switch SC
        case 0
            for tc = 1 : tx+100 
                time = tc * dt;
                tr = w * dt * tc;
                    if sweep == 1
                        freq = freq_start + freq_add*(time/cal_time);
                    else
                        freq = 10000;
                    end
                w = 2 * pai * freq ;
                j = sqrt(-1);
                if tc < pai / (w * dt)
                     pin(tc) = sin(w * time);
                else
                    if ippa == 1
                        pin(tc) = 0;
                    end
                    if ippa == 0
                        pin(tc) = sin(w * time);
                    end
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
                        freq = freq_start + freq_add*(time/cal_time);
                    else
                        freq = freq_abs;
                    end
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

    
% 時間ループ開始
if mode == 0
for t = 1: tx
   time = t * dt;
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

%   pinで定義した波動の印加と音圧更新式
    for i = 2 : ix + 1
        for j = 2 : jx + 1
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
                        elseif j >= jd_b && j <= jd_b2
                            p1(i,j) = pin(t + s);
                        elseif j >= jd_c && j <= jd_c2
                            p1(i,j) = pin(t + 2*s);
                        elseif j >= jd_d && j <= jd_d2
                            p1(i,j) = pin(t + 3*s);
                        elseif j >= jd_e && j <= jd_e2
                            p1(i,j) = pin(t + 4*s);
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

% 境界条件
    switch(boundary)
        case 0
            i = 1;
                for j = 2 : jx
                    p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
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
                end
            j = jx + 1;
                for i = 2 : ix
                    p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
                end
        case 1
            for i = 1 : 1
                for j = 1 : jx + 1
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            end
            i = ix + 1;
                for j = 1 : jx + 1
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            for j = 1 : 1
                for i = 1 : ix + 1
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
            end
            j = jx + 1 ;
                for i = 1 : ix + 1
                    u1(i,j) = 0;
                    v1(i,j) = 0;
                end
        case 2
            i = 1;
                for j = 2 : jx 
                    p1(i,j) = a .* (p1(i + 1,j)-p2(i,j)) - b .* (p1(i + 2,j) - 2 .* p2(i + 1,j) + p3(i,j)) - c .* (p2(i + 2,j) - p3(i + 1,j)) + d .* p2(i + 1,j) - e .* p3(i + 2,j);
                end
            i = ix + 1;
                for j = 1 : jx
                    p1(i,j) = a .* (p1(i - 1,j) - p2(i,j)) - b .* (p1(i-2,j) - 2*p2(i-1,j) + p3(i,j)) - c .* (p2(i - 2,j) - p3(i-1,j)) + d .* p2(i -1,j) - e .* p3(i-2,j);
                end
            j = 1;
                for i = 1 : ix
                    p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
                end
            j = jx + 1;
                for i = 2 : ix 
                    p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
                end
    end
    
%   音圧の時間的更新(p1>p2>p3の順で1ステップ新しく、次の音波計算のため時間を更新する)
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            p3(i,j) = p2(i,j);
            p2(i,j) = p1(i,j);
        end
    end
    
%   音圧更新
    for i = 1 : ix
        for j = 1 : jx + 1
            u1(i,j) = cv1 .* u2(i,j) - cv2 .* (p2(i + 1, j) - p2(i,j));
        end
    end
    for i = 1 : ix + 1
        for j = 1 : jx
            v1(i,j) = cv1 * v2(i,j) - cv2 * (p2(i, j + 1) - p2(i,j));
        end
    end
%   凹みによる完全反射
    if hekomi == 1 || hekomi == 2
        for i = b_x : b_x + h_x
            for j = 1 : w_x
                u1(i,j) = 0;
                v1(i,j) = 0;
            end
        end
    end
%   凹みが2点の場合
    if hekomi == 2
        for i = b_x_2 : b_x_2 + h_x_2
            for j = 1 : w_x_2
                u1(i,j) = 0;
                v1(i,j) = 0;
            end
        end
    end
    
%   階段近似を用いた境界条件
%   幾何をy=ax+b等の方程式に落とし込み実装
    if stair == 1
        grad = yrange / xrange;
        for x1_s = 0 : dx : xrange / 2
            y1_s = yrange/2 - grad*x1_s;
            i1 = round(x1_s/dx) + 1;
            j1 = round(y1_s/dx) + 1;
            u1(i1, j1) = 0;
            v1(i1, j1) = 0;
            if j1 < 3
                for jstair = 1 : j1 + 1
                    u1(i1, jstair) = 0;
                    v1(i1, jstair) = 0;
                end
            else
                for jstair = j1 - 2 : j1 + 1
                    u1(i1, jstair) = 0;
                    v1(i1, jstair) = 0;
                end
            end
        end
        for x1_s = xrange / 2 : dx : xrange
            y1_s = grad*x1_s - yrange/2;
            i1 = round(x1_s/dx) + 1;
            j1 = round(y1_s/dx) + 1;
            u1(i1, j1) = 0;
            v1(i1, j1) = 0; 
            if j1 < 3
                for jstair = 1 : j1 + 1
                    u1(i1, jstair) = 0;
                    v1(i1, jstair) = 0;
                end
            else
                for jstair = j1 - 2 : j1 + 1
                    u1(i1, jstair) = 0;
                    v1(i1, jstair) = 0;
                end
            end
        end
        for x1_s = xrange/2 : dx : xrange
            y1_s = 1.5*yrange - grad*x1_s;
            i1 = round(x1_s/dx) + 1;
            j1 = round(y1_s/dx) + 1;
            u1(i1, j1) = 0;
            v1(i1, j1) = 0;
            if j1 > jx - 2
                for jstair = j1 : jx + 1 
                    u1(i1, jstair) = 0;
                    v1(i1, jstair) = 0;
                end
            else
                for jstair = j1 : j1 + 3
                    u1(i1, jstair) = 0;
                    v1(i1, jstair) = 0;
                end
            end
        end
        for x1_s = 0 : dx : xrange / 2
            y1_s = grad*x1_s + yrange/2;
            i1 = round(x1_s/dx) + 1;
            j1 = round(y1_s/dx) + 1;
            u1(i1, j1) = 0;
            v1(i1, j1) = 0;
            if j1 > jx - 2
                for jstair = j1 : jx + 1 
                    u1(i1, jstair) = 0;
                    v1(i1, jstair) = 0;
                end
            else
                for jstair = j1 : j1 + 3
                    u1(i1, jstair) = 0;
                    v1(i1, jstair) = 0;
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
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
        end
        for istair = 3 : ix - 1
            for jstair = 1 : 2
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
            for jstair = jx  : jx + 1
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
        end
        for istair = ix : ix + 1
            for jstair = 1 : jx + 1
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
        end
    end
    if stair == 4
        for istair = 1 : 2
            for jstair = 1 : jx + 1
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
        end
        for istair = 3 : ix 
            nodo_length = yrange / 2;
            nodo_length = 0;
            kuchi_length = yrange;
            y_length = yrange;
            grad = (kuchi_length - nodo_length)/(2*(ix - 3)*dx);
            for jstair = 1 : round((0.5*y_length - nodo_length/2 - dx*grad*(istair - 3))/dx)
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
            for jstair = round((0.5*y_length + nodo_length/2 + dx*grad*(istair - 3))/dx) : jx + 1
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
        end
        for istair = ix : ix + 1
            for jstair = 1 : jx + 1
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
        end
    end
%   ロボットハンド
    if stair == 5
        theta = 40;
        div = tand(theta);
        j_0 = round(0.04/dx);
        j_1 = round((0.04 + 0.01*sind(theta))/dx);
        j_2 = round((0.04 + 0.07*cosd(theta))/dx);
        for istair = 1 : 2
            for jstair = 1 : jx + 1
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
        end
        for istair = 3 : ix - 1
            for jstair = 1 : 2
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
            for jstair = jx  : jx + 1
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
        end
        for istair = ix : ix + 1
            for jstair = 1 : jx + 1
                u1(istair, jstair) = 0;
                v1(istair, jstair) = 0;
            end
        end
%       1
        for j = 3 : round((0.04)/dx)
%             for i = 3 : round((0.01/dx))
%                 p1(i,j) = 0;
%             end
            for i = 1 : 2
                u1(i, j) = 0;
                v1(i, j) = 0;
            end
            for i = round((0.01/dx)) + 1 : round((0.01/dx)) + 8
                u1(i, j) = 0;
                v1(i, j) = 0;
            end
        end
%        2
        x1 = 3*dx;
        x2 = 0.01 - 0.01*cosd(theta);
        x3 = 0.01;
        x4 = x2 + 0.01/cosd(theta);
        for j = round(0.04/dx) : j_1
            xj = x1 + (j - j_0)*((x2 - x1)/(j_1 - j_0));
            xj_2 = x3 + (j - j_0)*((x4 - x3)/(j_1 - j_0));
            for i = round(xj/dx) - 2 : round(xj/dx) - 1
                u1(i, j) = 0;
                v1(i, j) = 0;
            end
            for i = round(xj_2/dx) + 1 : round(xj_2/dx) + 4
                u1(i, j) = 0;
                v1(i, j) = 0;
            end
        end
%        3
        x1 = 0.01 - 0.01*cosd(theta);
        x2 = (0.01 + 0.07*sind(theta)) - 0.01/cosd(theta);
        k = 0;
        for j = j_1 : j_2
            xj = x1 + (j - j_1)*((x2 - x1)/(j_2 - j_1));
%             for i = round(xj/dx) : round((xj + 0.01/cosd(theta))/dx)
%                 p1(i, j) = 0;
%             end
            for i = round(xj/dx) - 4 : round(xj/dx) - 1
                u1(i, j) = 0;
                v1(i, j) = 0;
            end
            for i = round((xj + 0.01/cosd(theta))/dx) + 1 : round((xj + 0.01/cosd(theta))/dx) + 4
                u1(i, j) = 0;
                v1(i, j) = 0;
            end
        end
    end
    
%   音圧の時間的更新（pと一緒）
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            u2(i,j) = u1(i,j);
            v2(i,j) = v1(i,j);
        end
    end
    
%   音圧の絶対値、ここでやる必要ないかも...笑
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            pressure(i,j) = (real(p1(i,j)) .^2).^0.5 ;
        end
    end
    
%   特定の位置の時間変化を配列に入れる
    if mod(t,5)==0
        p_keisoku_taihi(t/5) = p1(6, y_half);
%        p_keisoku_taihi(t/5) = p1(sokuteiten_x_g, sokuteiten_y_g);
    end
%   上のより詳細な情報を
    if mod(t, 1) == 0
        p_keisoku_detail(t) = p1(6, y_half);
    end
    
    
%   計算してる時間をコンソールに表示するとシミュレーション計算中安心して見れる
    if mod(t,500) == 0
        disp(t);
        time
    end
    
%   以降図へのプロット処理
    y = 1 : jx + 1 ;
    x = 1 : ix + 1;
    x_p = x * dx;
    y_p = y * dx;
    
%   位置-音圧プロット
    if mode_plot == 0
        if mod(t,100) == 0
            plot_image = pressure';
            imagesc(x_p, y_p, plot_image(y, x));
            %surf(x_p, y_p, abs(p1(x, y)));
            colorbar
            %colormap gray ;
%             title(['pressure when ',num2str(time),'seconds have passed'])
            title(['pressure when frequency =',num2str(freq),'Hz&&', num2str(time), '(s) have passed'])
            xlabel('x(mm)')
            ylabel('y(mm)')
            grid on;
            drawnow

            %アニメでなく画像で観たい時（savefig等にするとなおよし）        
            if time > snap_time1
                if snap1 == 0;
                    snap1 = snap1 + 1;
                    pause(2)
                end
            end
            if time > snap_time2
                if snap2 == 0;
                    snap2 = snap2 + 1;
                    pause(2)
                end
            end
            if time > snap_time3
                if snap3 == 0;
                    snap3 = snap3 + 1;
                    pause(2)
                end
            end
            if time > snap_time4
                if snap4 == 0;
                    snap4 = snap4 + 1;
                    pause(2)
                end
            end
        end
    end
    
%   x軸方向のみの位置-音圧プロット
    if mode_plot == 1
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
    
%   ホーン用（ほぼ使わない）
    if mode_plot ==2
        if time > 0.003
            for i = 1 : ix + 1
                px(i) = p1(i, y_half);
            end
            i = 1 : ix + 1;
            x_i = i * dx;
            plot(x_i, px)
            px = px';
            dlmwrite('nothorn_x.csv', px, 'precision', '%.10f', 'delimiter', ',')
            break;
        end
    end
    
%   時間-音圧プロット
    if mode_plot == 3
        if mod(t, 50) == 0
            t_x = 1 : t / 5;
            time = t_x *5 * dt;
            plot(time, p_keisoku_taihi(t_x));
            drawnow
        end
        if time > 0.01
                dlmwrite('kyokaisei.csv', p_keisoku_taihi, 'precision', '%.10f', 'delimiter', ',')
                break;
        end
    end
    
%   時間-音圧プロット
    if mode_plot == 4
        p_kei(t) = p1(1, y_half);
        p_kei(t)
        if time == cal_time - 0.001
            t_x = 1 : tx;
            plot(t_x * dt , p_kei(t_x))
            disp("計算が終了しました")
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
                p_keisoku_spec_col = p_keisoku_spec.';
                p_keisoku_taihi = p_keisoku_taihi.';
                p_keisoku_detail = p_keisoku_detail.';
                dlmwrite('pow600_1mm_0to12kHz_60ms_fin.csv', p_keisoku_taihi, 'precision', '%.15f', 'delimiter', ',')
                dlmwrite('pow600_1mm_0to12kHz_60ms_detail_fin.csv', p_keisoku_detail, 'precision', '%.15f', 'delimiter', ',')
                %dlmwrite('pow500_0to4_1cm.csv', p_keisoku_taihi, 'precision', '%.10f', 'delimiter', ',')
                break;
            end
    end
    
%   ある時間のx-pプロット
    if mode_plot == 7
        if time > 0.001
            plot(x_p,p1(x, y_half));
            break;
        end
    end
    
%   ある時間でのx-pプロット
    if mode_plot == 8
        sokutei_point = 300;
        if abs(p_keisoku_taihi(sokutei_point)) > 0
            sokutei_point * dx
            disp(sokutei_point * dx / time);
            break;
        end
        if mod(t,5) == 0
            plot(x_p,p1(x, y_half));
            grid on;
            drawnow
        end
    end
end
end

% 初期印加位置がどこかの検証（mode = 1以降時間ループは起こらず）
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

% 主要パラメータの表示
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

% 指向性シミュレータ
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

% 理論式のみのプロット　
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
    title(['pressure when frequency =',num2str(freq),'Hz&&', num2str(time), '(s) have passed'])
    xlabel('y(mm)')
    ylabel('x(mm)')
    grid on;
end

% 凹みを含んだ境界
if mode == 6
    for i = 1 : ix + 1
        for j = 1 : jx + 1
            p1(i, j) = 5;
        end
    end
        for istair = 1 : 2
            for jstair = 1 : jx + 1
                p1(istair, jstair) = 10;
            end
        end
        for istair = 3 : ix - 1
            for jstair = 1 : 2
                p1(istair, jstair) = 10;
            end
            for jstair = jx  : jx + 1
                p1(istair, jstair) = 10;
            end
        end
        for istair = ix : ix + 1
            for jstair = 1 : jx + 1
                p1(istair, jstair) = 10;
            end
        end
        p1(ix+1, jx+1) = 0;
        b_po = 0.3 ; %凹み位置
        h = 0.02;%凹み幅
        w = 0.001;%凹みふかさ
        b_x = round(b_po / dx) + 1;
        h_x = round(h / dx) + 1;
        w_x = round(w / dx) + 2;
        for i = b_x : b_x + h_x
            for j = 1 : w_x
                p1(i,j) = 10;
            end
        end
        x_i = 1 : ix + 1 ;
        y_i = 1 : jx + 1;
        x = x_i*dx;
        y = y_i*dx;
        p1 = p1';
        imagesc(x, y, p1(y_i,x_i));
        colorbar;
end

% ロボットハンドの初期印加位置と境界の位置表示
if mode == 7
    if range == 5
        theta = 40;
        div = tand(theta);
        j_0 = round(0.04/dx);
        j_1 = round((0.04 + 0.01*sind(theta))/dx);
        j_2 = round((0.04 + 0.07*cosd(theta))/dx);
        for i = 1 : ix + 1
            for j = 1 : jx + 1
                p1(i, j) = 5;
            end
        end
        for istair = 1 : 2
            for jstair = 1 : jx + 1
                p1(istair, jstair) = 10;
            end
        end
        for istair = 3 : ix - 1
            for jstair = 1 : 2
                p1(istair, jstair) = 10;
            end
            for jstair = jx  : jx + 1
                p1(istair, jstair) = 10;
            end
        end
        for istair = ix : ix + 1
            for jstair = 1 : jx + 1
                p1(istair, jstair) = 10;
            end
        end
%       1
        for j = 3 : round((0.04)/dx)
            for i = 3 : round((0.01/dx))
                p1(i,j) = 0;
            end
        end
%        2
        x1 = 3*dx;
        x2 = 0.01 - 0.01*cosd(theta);
        x3 = 0.01;
        x4 = x2 + 0.01/cosd(theta);
        for j = round(0.04/dx) : j_1
            xj = x1 + (j - j_0)*((x2 - x1)/(j_1 - j_0));
            xj_2 = x3 + (j - j_0)*((x4 - x3)/(j_1 - j_0));
            for i = round(xj/dx) : round(xj_2/dx)
                p1(i,j) = 0;
            end
        end
%        3
        x1 = 0.01 - 0.01*cosd(theta);
        x2 = (0.01 + 0.07*sind(theta)) - 0.01/cosd(theta);
        k = 0;
        for j = j_1 : j_2
            xj = x1 + (j - j_1)*((x2 - x1)/(j_2 - j_1));
            for i = round(xj/dx) : round((xj + 0.01/cosd(theta))/dx)
                p1(i, j) = 0;
            end
        end
        x_i = 1 : ix + 1 ;
        y_i = 1 : jx + 1;
        x = x_i*dx;
        y = y_i*dx;
        p1 = p1';
        imagesc(x, y, p1(y_i,x_i));
        colorbar;
    end
end
