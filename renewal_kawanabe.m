%%
% オシロスコープから抽出したデータを二回フーリエ変換して位置応答に変換する
%%%%オーバーラップあり

%% 設定

clear all;
close all;

set(0, 'DefaultUicontrolFontSize', 15);
set(0, 'defaultAxesFontName', 'メイリオ');
set(0, 'defaultTextFontName', 'メイリオ');
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 15);
set(0, 'DefaultAxesLineWidth', 1);
set(0, 'DefaultFigureColor', 'w');


%% 各種パラメーター

e = 0.1; % 入力電圧
opol = 100; % トレンド近似の字数を決定

m = 25;% 窓間隔
d = 1; % 分割間隔

f1 = 800; % スイープ開始周波数（Hz）
f2 = 7800; % 終了周波数

st = 1300; % フーリエ変換の開始周波数（Hz）
ed = 7300; % 終了周波数

load_data = csvread('hekoari06004000col.csv'); % 2行目より下を読み込む
noload_data = csvread('hekonashi1500col.csv'); % 2行目より下を読み込む


st_wave = 104;  % スイープ波形の始まる時間　可変
c = 340; % 大気中での音速 (m/s)
%% スイープ信号の表示



t = load_data(:,1); % 時間軸
v_1 = load_data(:,4)./e; % 入力電圧で割って定数化
v_2 = load_data(:,4)./e; % 入力電圧で割って定数化
v_o = noload_data(:,4)./e;


figure('Name', 'スイープ応答信号', 'NumberTitle', 'off')
plot(t*10^3, v_1, 'b');
xlabel('Time (ms)');
ylabel('Relative response (arb)');


%% 時間軸から周波数スペクトルを得る（負荷有りと負荷無し）


n_fft = 2^12; % 4096点FFT
dt = t(2) - t(1); % サンプリング周期
Fs = 1/dt; % 周波数領域での周期



st_1 = st_wave:d:st_wave+1800; % スイープ波形の分割
ed_1 = st_wave+m:d:st_wave+1800+m;


for i = 1:1:length(st_1); % オーバーラップ法にて周波数特性を作成
    
    st_2 = round(st_1(i)*10^-4/dt); % スイープ波形から切り出し開始
    ed_2 = round(ed_1(i)*10^-4/dt); % 切り出し終わり

    v_2 = v_1(st_2:ed_2); 
    v_3 = v_o(st_2:ed_2);

    Vl_f = (abs(fft(v_2, n_fft))).^2; % 切り出した波形をフーリエ変換
    Vo_f = (abs(fft(v_3, n_fft))).^2;

    peak_l(i) = max(Vl_f); % 最大値を記録していく
    peak_o(i) = max(Vo_f);

end 


f3 = f1:(f2 - f1)/(length(st_1)-1):f2; % 周波数軸の作成

figure('Name', '負荷有りの周波数特性', 'NumberTitle', 'off')
plot(f3*10^-3,peak_l,'b');
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');


figure('Name', '負荷無しの周波数特性', 'NumberTitle', 'off')
plot(f3*10^-3, peak_o, 'b');
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');


%% 移動平均を用いて平滑化処理（負荷有りと負荷無し）

    
for k = 0:1:4; % 移動平均で周波数特性を平滑化
    
    number = 4; % 各データ点付近の平均取る個数（偶数）

    zero = zeros(1, number/2); % ゼロ行列の作成

    peak_C = horzcat(zero, peak_l, zero); % 周波数特性の前後にゼロ行列を追加する
    peak_D = horzcat(zero, peak_o, zero);
    
for i = 1:1:length(st_1);
    
    peak_i = peak_C(i:i+number); % 周波数特性から切り取る
    peak_ave(i) = mean(peak_i); % 切り取った波形の平均を記録していく
    peak_l(i) = peak_ave(i); % 平均値を新しい行列に入れていく
    
    
    peak_i = peak_D(i:i+number);
    peak_ave(i) = mean(peak_i);
    peak_o(i) = peak_ave(i);
    
end

end


figure('Name', '平滑化処理後の周波数特性', 'NumberTitle', 'off')
plot(f3*10^-3, peak_o, 'r', 'linewidth', 1.5);

hold on

plot(f3*10^-3, peak_l, 'b', 'linewidth', 1.5);
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');


%% リプルを抽出しトレンドの除去


peak_l = peak_l - peak_o; % 負荷有りと負荷無しの差分を抽出

figure('Name', '直流成分除去前のリプル', 'NumberTitle', 'off')
plot(f3*10^-3, peak_l, 'b', 'linewidth', 2);
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');
%ylim([-3000 3000]);


[p, s, mu] = polyfit(f3, peak_l, opol); % リプルのトレンドを作成
f_y = polyval(p, f3, [], mu);

peak_l = peak_l - f_y; % トレンドを除去



figure('Name', '直流成分除去後のリプル', 'NumberTitle', 'off')
plot(f3*10^-3, peak_l, 'b', 'linewidth', 2);
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');
%ylim([-3000 3000]);


%% 周波数スペクトルから位置応答を検出


n = 2^14;
df3 = (f2 - f1)/(length(st_1) - 1); % 周波数領域のサンプリング間隔
dx = (c/2)/df3; % 位置領域での周期

st_3 = round((st - f1)/df3);
ed_3 = round((ed - f1)/df3);


ripple_f = peak_l(st_3:ed_3); % 波形の切り出し
% ripple_f = ripple_f.*hamming(length(ed_3 - st_3));



V_x = fft(ripple_f, n); % フーリエ変換して位置応答を作成


x = 0:dx/n:dx - dx/n; % 位置軸を作成
x = x';


for k = 0:1:2; % 位置応答を平滑化
    
    number = 4; % 各データ点付近の平均取る個数（偶数）

    zero = zeros(1, number/2);

    ave_response = horzcat(zero, V_x, zero);

    
for i = 1:1:length(V_x);
    
    V_x_i = ave_response(i:i+number);
    V_show(i) = mean(V_x_i);
    
    
end

end

noise = ones(n, 1).*0.8;

figure('Name', '位置応答', 'NumberTitle', 'off')
plot(x*10^3, abs(V_show)*10^-4, 'r', 'linewidth', 2);
xlabel('Position (mm)');
ylabel('Relative response (arb)');
xlim([0 2500]);



%% 応答ピークが負荷の影響かを判定


[peak_val, nx] = max(abs(V_x).*10^-4); % 位置応答のピークを抽出

peak_x = (nx - 1)*dx/n; % ピーク位置の場所を算出
peak_x_for_display = round(peak_x*10^3);


if peak_val > 3; % 閾値を用いてピークが負荷によるものかを判定
    display('負荷あり');
    display(peak_val);
    display(round(peak_x*10^3));
    
else display('負荷無し');
    display(peak_val);
    display(round(peak_x*10^3));
    
    
end
