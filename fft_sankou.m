%%
% �I�V���X�R�[�v���璊�o�����f�[�^����t�[���G�ϊ����Ĉʒu�����ɕϊ�����
%%%%�I�[�o�[���b�v����

clear all ;

c = 340; %����
c0 = 340;
% rou0 = 1.293; %���x�ikg/m^3
cal_time = 0.18;
rou0 = 1.293;
% rou0 = 1000;?
freq_param = 0.129/0.17;
freq_a = 2000;
freq_start = 1000*freq_param;
freq_add = 2000*freq_param;
freq = freq_a*freq_param;
freq_abs = freq_a*freq_param;
ramuda = c0 / freq;
dx_param = 0.01; %0.05-0.025 
dx = ramuda*dx_param; % �ɂ�20-30���̈�
% dt = dx / (5 * c);
% dt = dx*0.15 / (c0);
crn_param = 0.2;
dt = dx*crn_param/ (c0);
dt = 5*dt;
% �N�[����������𗧂Ă�]

%% �ݒ�

% clear all;
% close all;

set(0, 'DefaultUicontrolFontSize', 15);
set(0, 'defaultAxesFontName', '���C���I');
set(0, 'defaultTextFontName', '���C���I');
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 15);
set(0, 'DefaultAxesLineWidth', 1);
set(0, 'DefaultFigureColor', 'w');


%% �e��p�����[�^�[

e = 1; % ���͓d��
opol = 100; % �g�����h�ߎ��̎���������

m = 25;% ���Ԋu
d = 1; % �����Ԋu

% f1 = 800; % �X�C�[�v�J�n���g���iHz�j
% f2 = 7800; % �I�����g��
f1 = 1000;
f2 = 3000;

st = 1000; % �t�[���G�ϊ��̊J�n���g���iHz�j
ed = 3000; % �I�����g��

csvrangemax = cal_time/(5*dt) - mod(cal_time/(5*dt), 100);

load_data = csvread('kairyouheko0300.csv'); % 2�s�ڂ�艺��ǂݍ���
noload_data = csvread('kairyouhekonashi.csv'); % 2�s�ڂ�艺��ǂݍ���

load_data = load_data(1:csvrangemax);
noload_data = noload_data(1:csvrangemax);

st_wave = 1;  % �X�C�[�v�g�`�̎n�܂鎞�ԁ@��
c = 340; % ��C���ł̉��� (m/s)
%% �X�C�[�v�M���̕\��

cal_time_max_int = cal_time/(5*dt) - mod(cal_time/(5*dt), 100);
cal_time_max = cal_time_max_int * 5 * dt;
t = 5*dt : 5*dt : cal_time_max;

% t = load_data(:,1); % ���Ԏ�

% v_1 = load_data(:,4)./e; % ���͓d���Ŋ����Ē萔��
v_1 = load_data./e;
% v_2 = load_data(:,4)./e; % ���͓d���Ŋ����Ē萔��
v_2 = load_data./e;
% v_o = noload_data(:,4)./e;
v_o = noload_data./e;


figure('Name', '�X�C�[�v�����M��', 'NumberTitle', 'off')
plot(t*10^3, v_1, 'b');
xlabel('Time (ms)');
ylabel('Relative response (arb)');


%% ���Ԏ�������g���X�y�N�g���𓾂�i���חL��ƕ��ז����j


n_fft = 2^12; % 4096�_FFT
dt = t(2) - t(1); % �T���v�����O����
Fs = 1/dt; % ���g���̈�ł̎���



% st_1 = st_wave:d:st_wave+1800; % �X�C�[�v�g�`�̕���
% ed_1 = st_wave+m:d:st_wave+1800+m;
st_1 = st_wave:d:st_wave+1650; % �X�C�[�v�g�`�̕���
ed_1 = st_wave+m:d:st_wave+1650+m;


for i = 1:1:length(st_1); % �I�[�o�[���b�v�@�ɂĎ��g���������쐬
    
%     st_1(i)
    i;
    st_2 = round(st_1(i)*10^-4/dt); % �X�C�[�v�g�`����؂�o���J�n
    ed_2 = round(ed_1(i)*10^-4/dt); % �؂�o���I���
    
    v_2 = v_1(st_2:ed_2); 
    v_3 = v_o(st_2:ed_2);

    Vl_f = (abs(fft(v_2, n_fft))).^2; % �؂�o�����g�`���t�[���G�ϊ�
    Vo_f = (abs(fft(v_3, n_fft))).^2;

    peak_l(i) = max(Vl_f); % �ő�l���L�^���Ă���
    peak_o(i) = max(Vo_f);

end 


f3 = f1:(f2 - f1)/(length(st_1)-1):f2; % ���g�����̍쐬

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


[p, s, mu] = polyfit(f3, peak_l, opol); % ���v���̃g�����h���쐬
f_y = polyval(p, f3, [], mu);

peak_l_removed = peak_l - f_y; % �g�����h������



figure('Name', '��������������̃��v��', 'NumberTitle', 'off')
plot(f3*10^-3, peak_l_removed, 'b', 'linewidth', 2);
xlabel('Frequency (kHz)');
xlim([f1*10^-3 f2*10^-3]);
ylabel('Relative response (arb)');
%ylim([-3000 3000]);


%% ���g���X�y�N�g������ʒu���������o


n = 2^14;
df3 = (f2 - f1)/(length(st_1) - 1); % ���g���̈�̃T���v�����O�Ԋu
dx = (c/2)/df3; % �ʒu�̈�ł̎���

st_3 = round((st - f1)/df3);
ed_3 = round((ed - f1)/df3);

if st_3 == 0
    st_3 = 1;
end


ripple_f = peak_l(st_3:ed_3); % �g�`�̐؂�o��
% ripple_f = ripple_f.*hamming(length(ed_3 - st_3));



V_x = fft(ripple_f, n); % �t�[���G�ϊ����Ĉʒu�������쐬


x = 0:dx/n:dx - dx/n; % �ʒu�����쐬
x = x';


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
xlim([0 2500]);



%% �����s�[�N�����ׂ̉e�����𔻒�


[peak_val, nx] = max(abs(V_x).*10^-4); % �ʒu�����̃s�[�N�𒊏o

peak_x = (nx - 1)*dx/n; % �s�[�N�ʒu�̏ꏊ���Z�o
peak_x_for_display = round(peak_x*10^3);


if peak_val > 3; % 臒l��p���ăs�[�N�����ׂɂ����̂��𔻒�
    display('���ׂ���');
    display(peak_val);
    display(round(peak_x*10^3));
    
else display('���ז���');
    display(peak_val);
    display(round(peak_x*10^3));
    
    
end
