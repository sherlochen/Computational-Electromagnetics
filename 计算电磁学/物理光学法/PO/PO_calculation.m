clear all
close all
%% ***********************
% �������䲨
c0 = 3e8;%����
f = 10e9;%���䲨Ƶ��
lamda = c0/f;
K0 = 2*pi/(c0/f);%����
theta  = 15/180*pi;% ����Ƕ�����theta phi
phi = 10/180*pi;
Ki_x = sin(theta)*cos(phi);
Ki_y = sin(theta)*sin(phi);
Ki_z = -cos(theta);
Ki_n = [Ki_x,Ki_y,Ki_z]/norm([Ki_x,Ki_y,Ki_z]);
Ks_n = -Ki_n;
Ki_E = [1,0,0];
Ki = [Ki_x*K0,Ki_y*K0,Ki_z*K0];
%% ***********************
% ����ģ������
% v�������������εĶ���[3 * n x 3] 
% f��������ÿ����������Ķ����б�[n x 3] 
% n����ÿ����������ķ���[n x 3]
% c�ǿ�ѡ�ģ����Ұ���5λ[n x 3]�Ĳ�ɫrgb����
% stltitle����ָ����stl�ļ��ı���[1 x 80]
filename = "C:\Users\001\Desktop\�����ѧ��\����\sphere\wangge.stl";
% filename = "C:\Users\001\Desktop\�����ѧ��\����\cone\cone.stl";
% filename = "C:\Users\001\Desktop\ȺĿ��������һ�廯\cylinder\cylinder.stl";
[v, p, n] = stlread(filename);
figure
patch('Faces',p,'Vertices',v,'FaceColor','white')
view(3);
%% ***********************
% �ж���Ԫ�Ƿ��ڵ�
triangle_number = length(p); %��Ԫ����
triangle_O1 = zeros(triangle_number,1);  %��Ԫ�ж��ڵ����� 1Ϊ������0Ϊ����
triangle_C1 = zeros(3*triangle_number,1);
%***************************************
% ����һ�����ڵ����� �����䲨����ͬ��Ԫ���������
for i = 1:triangle_number
    triangle_c = n(i,1)*Ki_n(1,1) + n(i,2)*Ki_n(1,2) + n(i,3)*Ki_n(1,3);
    if(triangle_c>0||triangle_c==0)     %���䲨ͬ������Ԫ����������ǣ�����Ԫ��ֵΪ0���ڵ�
        triangle_O1(i)= 0;
        triangle_C1(1+(i-1)*3:3+(i-1)*3,:)=0;
    elseif(triangle_c<0) %���䲨ͬ������Ԫ�������ɶ۽ǣ�����Ԫ��ֵΪ1δ�ڵ�
        triangle_O1(i)= 1;
        triangle_C1(1+(i-1)*3:3+(i-1)*3,:)=1;
    end
end
% figure
% handle = patch('Faces',p,'Vertices',v,'FaceVertexCData',triangle_C1,'FaceColor','interp');
% set(handle,'EdgeColor','black','LineWidth',0.1);
% colorbar
% view(3);
%% ******************************
% �������Ԫ����������
triangle_Z_points = zeros(triangle_number,3);
triangle_R = zeros(triangle_number,1);
triangle_Rs = zeros(triangle_number,1);
triangle_S = zeros(triangle_number,1);
for i = 1:triangle_number
    triangle_L1 = sqrt((v(p(i,2),1)-v(p(i,3),1))^2+(v(p(i,2),2)-v(p(i,3),2))^2+(v(p(i,2),3)-v(p(i,3),3))^2);
    triangle_L2 = sqrt((v(p(i,1),1)-v(p(i,3),1))^2+(v(p(i,1),2)-v(p(i,3),2))^2+(v(p(i,1),3)-v(p(i,3),3))^2);
    triangle_L3 = sqrt((v(p(i,1),1)-v(p(i,2),1))^2+(v(p(i,1),2)-v(p(i,2),2))^2+(v(p(i,1),3)-v(p(i,2),3))^2);
    triangle_Z_points(i,1)=(triangle_L1*v(p(i,1),1)+triangle_L2*v(p(i,2),1)+triangle_L3*v(p(i,3),1))/(triangle_L1+triangle_L2+triangle_L3);
    triangle_Z_points(i,2)=(triangle_L1*v(p(i,1),2)+triangle_L2*v(p(i,2),2)+triangle_L3*v(p(i,3),2))/(triangle_L1+triangle_L2+triangle_L3);
    triangle_Z_points(i,3)=(triangle_L1*v(p(i,1),3)+triangle_L2*v(p(i,2),3)+triangle_L3*v(p(i,3),3))/(triangle_L1+triangle_L2+triangle_L3);
    triangle_R(i) =  dot(Ki_n,triangle_Z_points(i,:));   % �������Ԫ���ĵ�����ԭ��ľ��������䷽���ϵ�ͶӰ
    triangle_Rs(i) = dot(triangle_Z_points(i,:),Ks_n);  % �������Ԫ���ĵ�����ԭ��ľ�����ɢ�䷽���ϵ�ͶӰ
    triangle_P = (triangle_L1+triangle_L2+triangle_L3)/2;  
    triangle_S(i) = sqrt(triangle_P*(triangle_P-triangle_L1)*(triangle_P-triangle_L2)*(triangle_P-triangle_L3));
end
figure
handle = patch('Faces',p,'Vertices',v,'FaceVertexCData',triangle_C1,'FaceColor','interp');
set(handle,'EdgeColor','black','LineWidth',0.1);
% colorbar
colormap gray
view(3);
hold on
scatter3(triangle_Z_points(:,1),triangle_Z_points(:,2),triangle_Z_points(:,3),3,'black','filled');
hold off

%% ******************************
% ����Ӧ����
% �������䲨�ﵽ����Ԫ�ĵ糡���ų�
% A_E = 1;%����糡����
Z0 = 120*pi; %���ɿռ��еĲ��迹
triangle_E = zeros(triangle_number,3);  %�������Ԫ�ĵ糡
triangle_H = zeros(triangle_number,3);
for i = 1:triangle_number
    triangle_E(i,:) = exp(-1i*K0*triangle_R(i)).*Ki_E;
    triangle_H(i,:) = cross(Ki_n,triangle_E(i,:));
end
% ��������Ӧ���� Jn=2n*Hin, �ڰ�����Ӧ����Ϊ0
triangle_J = zeros(triangle_number,3);
triangle_Es = zeros(triangle_number,3);
for i = 1:triangle_number
    if(triangle_O1(i)==1)    
          triangle_J(i,:) = cross(n(i,:),triangle_H(i,:));
          triangle_Es(i,:) = cross(triangle_J(i,:),Ks_n)*exp(1i*K0*triangle_Rs(i))*triangle_S(i);
    end
end
R = 100e10;
Es_1 = sum(triangle_Es);
ES1 = -1i*K0*exp(-1i*K0*R)/(2*pi*R)*cross(Ks_n,Es_1);
RCS1 = 10*log10(4*pi*R^2*(abs(ES1(1,1))+abs(ES1(1,2))+abs(ES1(1,3)))^2);
%% ***************************************
triangle_Rcs = zeros(triangle_number,1);
for i = 1:triangle_number
    if(triangle_O1(i)==1)    
         triangle_Rcs(i) = dot(n(i,:),Ks_n)*exp(2*1i*K0*dot(Ks_n,triangle_Z_points(i,:)))*triangle_S(i);
    end
end
ES2 = -1i*K0*exp(-1i*K0*R)/(2*pi*R)*sum(triangle_Rcs).*Ki_E;
RCS0 = 10*log10(4*pi*(abs(ES2(1,1))+abs(ES2(1,2))+abs(ES2(1,3)))^2);
RCS2 = 10*log10(4*pi/(lamda^2)*(abs(sum(triangle_Rcs)))^2);
