%%y
% �I�V���X�R�[�v���璊�o�����f�[�^����t�[���G�ϊ����Ĉʒu�����ɕϊ�����
%%%%�I�[�o�[���b�v����
clear
bai = 1;
for i = 1 : 2

c = 340; %����
c0 = 340;
% rou0 = 1.293; %���x�ikg/m^3
long = 3;
if long == 1
    cal_time = 0.18;
end
if long == 0
    cal_time = 0.06;
end
if long == 2
    cal_time = 0.36;
end
if long == 3
    cal_time = 0.01;
end
rou0 = 1.293;
% rou0 = 1000;
freq_param = 1;
freq_a = 1000;
freq_start = 1000;
freq_add = 2000;
freq = freq_a;
freq_abs = freq_a*freq_param;
ramuda = c0 / freq;
dx_param = 0.005; %0.05-0.025 
dx = ramuda*dx_param; % �ɂ�20-30���̈�
dx = 0.0002;
% dt = dx / (5 * c);
% dt = dx*0.15 / (c0);
crn_param = 0.2;
dt = dx*crn_param/ (c0);
% �N�[����������𗧂Ă�]

%% �ݒ�


set(0, 'DefaultUicontrolFontSize', 15);
set(0, 'defaultAxesFontName', '���C���I');
set(0, 'defaultTextFontName', '���C���I');
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 15);
set(0, 'DefaultAxesLineWidth', 1);
set(0, 'DefaultFigureColor', 'w');

%% �e��p�����[�^�[

e = 0.1; % ���͓d��
opol = 1000; % �g�����h�ߎ��̎���������

m = 25;% ���Ԋu
d = 1; % �����Ԋu

xL = 0.3;
achieve_time = 2*xL/c0;
start_time_g = round(achieve_time / (5*dt));

f1 = 000; % �X�C�[�v�J�n���g���iHz�j
f2 = 12000; % �I�����g��

st = 2000; % �t�[���G�ϊ��̊J�n���g���iHz�j
ed = 10000; % �I�����g��
% csvrangemax = cal_time/(5*dt);
% csvrangemax = cal_time/(dt) - mod(cal_time/(dt), 100)

load_data = csvread('robot_6cm_0.2mm_0to12kHz_10ms_20deg_fin.csv'); % 2�s�ڂ�艺��ǂݍ���
noload_data = csvread('robot_noload_0.2mm_0to12kHz_10ms_20deg_fin.csv'); % 2�s�ڂ�艺��ǂݍ���
csvrangemax = round(cal_time/(5*dt));
load_data = real(load_data(1 :csvrangemax));
noload_data = real(noload_data(1:csvrangemax ));

data_sabun = load_data - noload_data;

dlmwrite('sabun900.csv', data_sabun, 'precision', '%.15f', 'delimiter', ',')

dt_p = 5*dt;
t_sec = dt_p : dt_p : cal_time;
plot_start = round(0.0/dt_p)+1;
plot_end = round(cal_time/dt_p);

plot(t_sec(plot_start:plot_end), data_sabun(plot_start:plot_end))

% st_wave = 104;  % �X�C�[�v�g�`�̎n�܂鎞�ԁ@��
st_wave = 1;
c = 340; % ��C���ł̉��� (m/s)
%% �X�C�[�v�M���̕\��

cal_time_max_int = cal_time/(dt) - mod(cal_time/(dt), 100);
cal_time_max = cal_time_max_int * dt;
t_sec = 5*dt : 5*dt : cal_time;

t = 5 * dt : 5 * dt : cal_time; % ���Ԏ�
v_1 = load_data./e; % ���͓d���Ŋ����Ē萔��
v_2 = load_data./e; % ���͓d���Ŋ����Ē萔��
v_o = noload_data./e;

figure('Name', '�X�C�[�v�����M��', 'NumberTitle', 'off')
plot(t*10^3, v_1, 'b');
disp(size(t_sec))
disp(size(v_1))
plot(t_sec*10^3, v_1, 'b');
xlabel('Time (ms)');
ylabel('Relative response (arb)');


%% ���Ԏ�������g���X�y�N�g���𓾂�i���חL��ƕ��ז����j

n_fft = 2^12; % 4096�_FFT
dt = t(2) - t(1); % �T���v�����O����
Fs = 1/dt; % ���g���̈�ł̎���

if long == 1
    st_1 = st_wave:d:st_wave+1700; % �X�C�[�v�g�`�̕���
    ed_1 = st_wave+m:d:st_wave+1700+m;
end
if long == 0
    st_1 = st_wave:d:st_wave+570; % �X�C�[�v�g�`�̕���
    ed_1 = st_wave+m:d:st_wave+570+m;
end
if long == 2
    st_1 = st_wave:d:st_wave+3550; % �X�C�[�v�g�`�̕���
    ed_1 = st_wave+m:d:st_wave+3550+m;
end
if long == 3
%     st_1 = st_wave:d:st_wave+10000; % �X�C�[�v�g�`�̕���
%     ed_1 = st_wave+m:d:st_wave+10000+m;
    st_1 = st_wave:d:st_wave+970; % �X�C�[�v�g�`�̕���
    ed_1 = st_wave+m:d:st_wave+970+m;
end


for i = 1:1:length(st_1); % �I�[�o�[���b�v�@�ɂĎ��g���������쐬
    
    st_2 = round(st_1(i)*10^-5/dt); % �X�C�[�v�g�`����؂�o���J�n
    ed_2 = round(ed_1(i)*10^-5/dt); % �؂�o���I���

    v_2 = v_1(st_2:ed_2); 
    v_3 = v_o(st_2:ed_2);

    Vl_f = (abs(fft(v_2, n_fft))).^2; % �؂�o�����g�`���t�[���G�ϊ�
    Vo_f = (abs(fft(v_3, n_fft))).^2;

    peak_l(i) = max(Vl_f); % �ő�l���L�^���Ă���
    peak_o(i) = max(Vo_f);

end 

f3 = f1:(f2 - f1)/(length(st_1)-1):f2; % ���g�����̍쐬
length(f3)
length(peak_l)

figure('Name', '���חL��̎��g������', 'NumberTitle', 'off')
plot(f3*10^-3,peak_l,'b');
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');

figure('Name', '���ז����̎��g������', 'NumberTitle', 'off')
plot(f3*10^-3, peak_o, 'b');
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');
%% �ړ����ς�p���ĕ����������i���חL��ƕ��ז����j

    
for k = 0:1:4; % �ړ����ςŎ��g�������𕽊���
    
    number = 4; % �e�f�[�^�_�t�߂̕��ώ����i�����j

    zero = zeros(1, number/2); % �[���s��̍쐬

    peak_C = horzcat(zero, peak_l, zero); % ���g�������̑O��Ƀ[���s���ǉ�����
    peak_D = horzcat(zero, peak_o, zero);
    
for i = 1:1:length(st_1);
    
    peak_i = peak_C(i:i+number); % ���g����������؂���
    peak_ave(i) = mean(peak_i); % �؂������g�`�̕��ς��L�^���Ă���
    peak_l(i) = peak_ave(i); % ���ϒl��V�����s��ɓ���Ă���
    
    peak_i = peak_D(i:i+number);
    peak_ave(i) = mean(peak_i);
    peak_o(i) = peak_ave(i);
    
end

end

figure('Name', '������������̎��g������', 'NumberTitle', 'off')
plot(f3*10^-3, peak_o, 'r', 'linewidth', 1.5);

hold on

plot(f3*10^-3, peak_l, 'b', 'linewidth', 1.5);
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');


%% ���v���𒊏o���g�����h�̏���

peak_l = peak_l - peak_o; % ���חL��ƕ��ז����̍����𒊏o

figure('Name', '�������������O�̃��v��', 'NumberTitle', 'off')
plot(f3*10^-3, peak_l, 'b', 'linewidth', 2);
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');
%ylim([-3000 3000]);

peak_taihi = peak_l ;
[p, s, mu] = polyfit(f3, peak_l, opol); % ���v���̃g�����h���쐬
f_y = polyval(p, f3, [], mu);


peak_l = peak_l - f_y; % �g�����h������


figure('Name', '��������������̃��v��', 'NumberTitle', 'off')
plot(f3*10^-3, peak_l, 'b', 'linewidth', 2);
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');
%ylim([-3000 3000]);


%% ���g���X�y�N�g������ʒu���������o


n = 2^14;
df3 = (f2 - f1)/(length(st_1) - 1); % ���g���̈�̃T���v�����O�Ԋu
dx = (c/2)/df3 ; % �ʒu�̈�ł̎���

st_3 = round((st - f1)/df3);
ed_3 = round((ed - f1)/df3);


ripple_f = peak_taihi(st_3:ed_3); % �g�`�̐؂�o��
ripple_f = ripple_f.*hamming(length(ed_3 - st_3));



V_x = fft(ripple_f, n); % �t�[���G�ϊ����Ĉʒu�������쐬
%V_x = fft(peak_l, n); % �t�[���G�ϊ����Ĉʒu�������쐬

x = 0:dx/n:dx - dx/n; % �ʒu�����쐬
x = x';

if bai == 1
    x = x / 2;
end


for k = 0:1:2; % �ʒu�����𕽊���
    
    number = 4; % �e�f�[�^�_�t�߂̕��ώ����i�����j

    zero = zeros(1, number/2);

    ave_response = horzcat(zero, V_x, zero);

    
for i = 1:1:length(V_x);
    
    V_x_i = ave_response(i:i+number);
    V_show(i) = mean(V_x_i);

end

end

noise = ones(n, 1).*0.8;

figure('Name', '�ʒu����', 'NumberTitle', 'off')

plot(x*10^3, abs(V_show)*10^-4, 'r', 'linewidth', 2);



xlabel('Position (mm)');
ylabel('Relative response (arb)');
xlim([0 130]);

%% �����s�[�N�����ׂ̉e�����𔻒�


[peak_val, nx] = max(abs(V_x).*10^-4); % �ʒu�����̃s�[�N�𒊏o

peak_x = (nx - 1)*dx/n; % �s�[�N�ʒu�̏ꏊ���Z�o
peak_x_for_display = round(peak_x*10^3);

if bai == 1
    peak_x = peak_x / 2;
end

if peak_val > 3; % 臒l��p���ăs�[�N�����ׂɂ����̂��𔻒�
    display('���ׂ���');
    display(peak_val);
    display(round(peak_x*10^3));
    
else display('���ז���');
    display(peak_val);
    display(round(peak_x*10^3));
end
break;
end
