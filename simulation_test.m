xrange = 2.02;
yrange = 0.02;
dx = 0.001;
dt =   0.000002;
cal_time = 0.01;
x1 = 0.03;
x2 = 0.04;
y1 = 0.1;
y2 = 0.44;
yd = 0.008;
yd2 = 0.012;
xd = 0.001;
xd2 = 0.002;
t1 = 0;
t2 = 0;
speed = 0;
disp_hensu = 0;
f = 1000; %���g��
c0 = 340; %����
rou0 = 1.293; %���x�ikg/m^3
absp0 = - 0.5; % �z���W��
gensui0 = (f*absp0) / (8.686*c0); % �����W��
ix = round(xrange / dx); %x��Ԋ��o�̐�
jx = round(yrange / dx); %y��Ԋ��o�̐�
tx = fix(cal_time / dt ); %���Ԋ��o�̐�
td = 15; % ���g���̑���i�H�j
id = round(xd / dx);%x�̈ʒu�i�H�j
id2 = round(xd2 / dx ) ;
jd = round(yd / dx) ; %y��Ԋ��o�̑���i�H�j
jd2 = round(yd2/dx);
snap = 10 ;%�����t�@�C���̏o�͊��o
del_T = round((cal_time/dt)/snap);%�X�i�b�v�V���b�g�̎��Ԋ��o�̑��
ix1 = round(x1/dx);
ix2 = round(x2/dx);
jx1 = round(y1/dx);
jx2 = round(y2/dx);
% p1 = zeros(ix+1,jx+1);
% p2 = zeros(ix+1,jx+1);
p3 = zeros(ix+1,jx+1);
u1 = zeros(ix+1,jx+1);
u2 = zeros(ix+1,jx+1);
v1 = zeros(ix+1,jx+1);
v2 = zeros(ix+1,jx+1);
mo = zeros(ix+1,jx+1); %�`���[�u���f��
pressure = zeros(ix + 1 , jx + 1);
SC = 0 ;%��U�֐� �O�Ȃ�A��1�Ȃ�K�E�V�A��2�n�j���O3�����g���g
WN = 1;
W_end = round(2*WN/(f*dt)) - 1 ;
pin = zeros(1,tx);
crn =(c0 * dt)/dx ; %�N�[������
% crn = 0.2;
dd = 199/200 ;%higdons absorption boundary
pai = 3.1415 ;
hasu_o0 = 2*pai*f/c0 ;%���ۂ̔g��
hasu0 = sqrt(hasu_o0*hasu_o0 - gensui0^2);%��������ꍇ�̔g��
c_m0 = 2*pai*f/hasu0;%�����̂���ꍇ�̉���
alpha0 = 2*hasu_o0*rou0*c_m0/hasu0;%�z����
kap0 = c_m0^2*rou0;
% cp1 = 1;
cp1 = 1;
% cp2 = kap0*dt/dx;
cp2 = rou0*c_m0^2*dt/dx;
% cv1 = (2*rou0 - alpha0*dt)/(2*rou0 + alpha0*dt);
cv1 = 1;
% cv2 = dt*2 /((2*rou0 + alpha0*dt)*dx);
cv2 = dt / (rou0 * dx);
a = 2*(crn - 1)/(crn + 1);
b = ((crn - 1)/(crn + 1))*2;%�v�����v����
c = 2 * (crn - 1) / (crn + 1) * dd;
d = dd * 2;
e = dd *2; %�v�����v����
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
%%�������@
for i = 1 : ix + 1
    for j = 1 : jx + 1
        p1(i,j) = 0;
        p2(i,j) = 0;
        p3(i,j) = 0;
        u1(i,j) = 0;
        u2(i,j) = 0;
        v1(i,j) = 0;
        v2(i,j) = 0;
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
w = 2 * pai * f ;
% switch SC
%     case 0
%         for t = 1 : tx
%             tr = w * dt * t / WN;
%             if tr <= pai
%                 pin(t) = sin((w * dt * t));
%             else
%                 pin(t) = sin(w * dt * t);
%             end
%         end
%     case 1
%         tau = WN / f;
%         for t = 1:W_end
%             al = (4/tau) * (4/tau);
%             pin =exp(-al * (dt*t - tau)*(dt*t - tau)) * sin(w*(dt*t - tau));
%         end
%     case 2
%         for t = 1 : W_end
%             pin = (1-cos(w*dt*t/WN)) * sin(w*dt*t)/2;
%         end
%     case 3
%         for t = 1 : W_end
%         pin = sin(w * dt * real(t - 1));
%         end
% end

p1 , p2 = p_in(SC,ix,jx,id,id2,jd,jd2,f,tx,dt,WN,t);
disp("�s��̃T�C�Y��");
disp(mean2(p1));
disp(mean2(p2));
disp(size(p1));
disp(size(p2));
disp(size(u2));
disp(size(v2));
%%timeloop
for t = 1: tx
    break;
   time = t * dt * 1.47;
   if crn >= 1
       disp("�N�[���������s�K�؂ł��B");
       break;
   end
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
    j = 1;
    for i = 2 : ix
        p1(i,j) = a .* (p1(i,j+1) - p2(i,j)) - b .* (p1(i,j+2) - 2 * p2(i,j+1) + p3(i,j)) - c .* (p2(i,j+2) - p3(i,j + 1)) + d .* p2(i,j + 1) - e .* p3(i,j + 2);
    end
    j = jx + 1;
    for i = 2 : ix
        p1(i,j) = a .* (p1(i,j-1) - p2(i,j)) - b .* (p1(i,j-2) - 2 * p2(i,j-1) + p3(i,j)) - c .* (p2(i,j-2) - p3(i,j - 1)) + d .* p2(i,j - 1) - e .* p3(i,j - 2);
    end
    
    p_in(SC,ix,jx,id,id2,jd,jd2)
    
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
    disp(t);
    y = 1 : jx + 1 ;
    x = 1 : ix + 1;
     imagesc(y,x,pressure(x,y));
        title(['pressure when ',num2str(time),'seconds have passed'])
        xlabel('y(mm)')
        ylabel('x(mm)')
     grid on;
     drawnow
end

function [p1, p2] = p_in(SC,ix,jx,id,id2,jd,jd2,f,tx,dt,WN,t)
    pin = p_in_f(SC, f, tx, dt, WN);
    for i = id : id2
        for j = jd : jd2
            p1(i,j) = pin(t);
            p2(i,j) = pin(t);
        end
    end
end

function pin = p_in_f(SC, f, tx, dt, WN)
    pai = 3.1415;
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
end
    
            
            
