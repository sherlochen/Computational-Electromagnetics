clear all
close all
% 参数定义
c0 = 3e8;%光速
f = 1e9;%入射波频率
lamda = c0/f;
K0 = 2*pi/(c0/f);%波数
theta  = 10;%入射角度设置theta phi
phi = 0;
Ki_x = sin(theta/180*pi)*cos(phi/180*pi);
Ki_y = sin(theta/180*pi)*sin(phi/180*pi);
Ki_z = cos(theta/180*pi);
Ki_n = [Ki_x,Ki_y,Ki_z]/norm([Ki_x,Ki_y,Ki_z]);%入射方向的单位矢量
Ks_n = -Ki_n;%后向散射方向矢量
Z0 = 120*pi; %自由空间中的波阻抗