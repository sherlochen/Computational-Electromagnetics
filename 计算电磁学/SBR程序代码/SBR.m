clear all
% close all
filename = "F:\������ѧ\SBR�������\5qiu_wangge.stl";
%% ������Ԥ����
[v, p, n] = stlread(filename);
% ���䲨��������
theta = 90/180*pi;
phi = 0/180*pi;
Ki_x = sin(theta)*cos(phi);
Ki_y = sin(theta)*sin(phi);
Ki_z = cos(theta);
Ki_n = [Ki_x,Ki_y,Ki_z]/norm([Ki_x,Ki_y,Ki_z]);
Ks_n = -Ki_n;
% *************************************
triangle_number = length(p); %��Ԫ����
% �ж���Ԫ�Ƿ��ڵ�
triangle_O1 = zeros(triangle_number,1);  %��Ԫ�ж��ڵ����� 1Ϊ������0Ϊ����
triangle_C1 = zeros(3*triangle_number,1);
k1 = 1;
for i = 1:triangle_number
    % ********************************************************************
    % ���ڵ����� �����䲨����ͬ��Ԫ��������� �ж������Ͱ���
    triangle_c = n(i,1)*Ki_n(1,1) + n(i,2)*Ki_n(1,2) + n(i,3)*Ki_n(1,3);
    if(triangle_c>0)     %���䲨ͬ������Ԫ����������ǣ�����Ԫ��ֵΪ0���ڵ�
        triangle_O1(i,1)= 0;
        triangle_C1(1+(i-1)*3:3+(i-1)*3,1)=0;
    elseif(triangle_c<0) %���䲨ͬ������Ԫ�������ɶ۽ǣ�����Ԫ��ֵΪ1δ�ڵ�
        triangle_O1(i,1)= 1;
        bright_triangle(k1,:) = p(i,:);
        bright_triangle_n(k1,:) = n(i,:);
        triangle_C1(1+(i-1)*3:3+(i-1)*3,1)=1;
        k1 = k1+1;
    end
end
figure
handle = patch('Faces',p,'Vertices',v,'FaceVertexCData',triangle_C1,'FaceColor','interp');
set(handle,'EdgeColor','black','LineWidth',0.1);
colorbar
view(3);
% ********************************************
%% ��ͶӰ��
%ȷ�������������
Y_theta_projection = zeros(length(v),1);
Z_phi_projection = zeros(length(v),1);
for i = 1:length(v)
     Y_theta_projection(i) = cos(theta)*cos(phi)*v(i,1) + cos(theta)*sin(phi)*v(i,2) -sin(theta)*v(i,3);
     Z_phi_projection(i) = -sin(phi)*v(i,1) + cos(phi)*v(i,2);
end
Y_theta_max = max(Y_theta_projection);
Y_theta_min = min(Y_theta_projection);
Z_phi_max = max(Z_phi_projection);
Z_phi_min = min(Z_phi_projection);
% ***************************************
%% �����������������
c0 = 3e8;
f = 10e9;%linspace(1e9,2e9,11);
lamda = c0./f;
deta_lamda = lamda/10;
Z0 = 120*pi; %���ɿռ䲨�迹
K = (2*pi)/lamda;
M = round((Y_theta_max-Y_theta_min)/deta_lamda);
N = round((Z_phi_max-Z_phi_min)/deta_lamda);
Y_theta = linspace(Y_theta_min,Y_theta_max,M);
Z_phi = linspace(Z_phi_min,Z_phi_max,N);
Distance_r_max = sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2);
Distance_r = -1.5*max(Distance_r_max)*ones(M,N);
point_xyz = zeros(M*N,3);
F = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);-sin(phi) cos(phi) 0];
bright_triangle_number = length(bright_triangle);
normal_t = zeros(length(bright_triangle),1);
normal_u = zeros(length(bright_triangle),1);
normal_v = zeros(length(bright_triangle),1);
ray_edgepoint_plane_number = zeros(M,N);
% ***********************************************
%% �����ߺ���Ԫ�ཻ�жϡ�
point_cell{M,N} = zeros;
for i = 1:M
    for j = 1:N
%         % �����������и��ǵ㱣��
          point_xyz = (inv(F)*[Distance_r(i,j);Y_theta(i);Z_phi(j)]).' ;
%         % �жϸ��ǵ��ཻ����Ԫ
          point_cell{i,j} = point_cell_import(v,bright_triangle,point_xyz,Ki_n);
    end
end
% ****************************************
%% �����߹ܹ��ࡿ
% ������깹��һ�����߹ܣ�ǰ�ĸ�����Ϊ��ʱ��仯�����߹ܽǵ㣬���һ������Ϊ���ĵ�
Ray_tube = zeros(5*(M-1)*(N-1),7); 
Ray_tube_edgepoint = zeros(4*(M-1)*(N-1),7);
Ray_tube_midepoint = zeros((M-1)*(N-1),7);
k2  = 0;
for i = 1:M-1
    for j = 1:N-1
        Ray_tube(4*k2+1,:) = point_cell{i,j};
        Ray_tube(4*k2+2,:) = point_cell{i+1,j};
        Ray_tube(4*k2+3,:) = point_cell{i+1,j+1};
        Ray_tube(4*k2+4,:) = point_cell{i,j+1};
        Ray_tube(5*k2+5,1:3) = (point_cell{i,j}(1,1:3)+point_cell{i+1,j+1}(1,1:3))/2;   
        Ray_tube_Intermediate_variables(1,:) = point_cell_import(v,bright_triangle,Ray_tube(5*k2+5,1:3),Ki_n);
        Ray_tube(5*k2+5,4:7) = Ray_tube_Intermediate_variables(1,4:7);
        
        Ray_tube_edgepoint(4*k2+1,:) = point_cell{i,j};
        Ray_tube_edgepoint(4*k2+2,:) = point_cell{i+1,j};
        Ray_tube_edgepoint(4*k2+3,:) = point_cell{i+1,j+1};
        Ray_tube_edgepoint(4*k2+4,:) = point_cell{i,j+1};
        
        Ray_tube_midepoint(k2+1,:) = (point_cell{i,j}+point_cell{i+1,j+1})/2;
        k2 = k2 + 1;
    end
end
%********************************************************************
k3 = 1;
k4 = 1;
for i = 1:length(Ray_tube_edgepoint(:,1))
    if Ray_tube_edgepoint(i,4) == 0
       Ray_tube_edgepoint_1(k3,:) = Ray_tube_edgepoint(i,:);
       k3 = k3 +1;
    else
       Ray_tube_edgepoint_2(k4,:) = Ray_tube_edgepoint(i,:);
       k4 = k4 +1;
    end
end
% **********************************************
figure
handle = patch('Faces',p,'Vertices',v,'FaceVertexCData',triangle_C1,'FaceColor','interp');
set(handle,'EdgeColor','black','LineWidth',0.1);
colorbar
view(3);
hold on
scatter3(Ray_tube_edgepoint_1(:,1),Ray_tube_edgepoint_1(:,2),Ray_tube_edgepoint_1(:,3),3,'black','filled');
scatter3(Ray_tube_edgepoint_2(:,1),Ray_tube_edgepoint_2(:,2),Ray_tube_edgepoint_2(:,3),3,'red','filled');
hold off
% ************************************
%% ���ж����߹ܵ���Ч�ԡ�
% Ray_tube_usefull_1 �� Ray_tube ÿ����7λ����ɣ�ǰ����Ϊ�����λ�ã�����λΪ�ཻ��Ԫ������λΪ�����λ��
Ray_tube_number = length(Ray_tube(:,1))/5; %���߹�����
k5 = 1; %��Ч���߹�����
for j =1:Ray_tube_number
       w1 = Ray_tube(5*(j-1)+1,4); % �ǵ����߶�Ӧ�ŵ���Ԫ
       w2 = Ray_tube(5*(j-1)+2,4);
       w3 = Ray_tube(5*(j-1)+3,4);
       w4 = Ray_tube(5*(j-1)+4,4);
       w5 = Ray_tube(5*(j-1)+5,4);
       if((w1*w2*w3*w4*w5)~=0)
           w12 = dot(bright_triangle_n(w1,:),bright_triangle_n(w2,:));
           w13 = dot(bright_triangle_n(w1,:),bright_triangle_n(w3,:));
           w14 = dot(bright_triangle_n(w1,:),bright_triangle_n(w4,:));
           w23 = dot(bright_triangle_n(w2,:),bright_triangle_n(w3,:));
           w24 = dot(bright_triangle_n(w2,:),bright_triangle_n(w4,:));
           w34 = dot(bright_triangle_n(w3,:),bright_triangle_n(w4,:));
           if(w12>=0.8&&w13>=0.8&&w14>=0.8&&w23>=0.8&&w23>=0.8&&w34>=0.8)
               % �������߹�׷����Ȳ�һ�£���ɵ��������߹��ݲ�����
               Ray_tube_usefull_1(5*(k5-1)+1,:) = Ray_tube(5*(j-1)+1,:); %��Ч���߹�
               Ray_tube_usefull_1(5*(k5-1)+2,:) = Ray_tube(5*(j-1)+2,:);
               Ray_tube_usefull_1(5*(k5-1)+3,:) = Ray_tube(5*(j-1)+3,:);
               Ray_tube_usefull_1(5*(k5-1)+4,:) = Ray_tube(5*(j-1)+4,:);
               Ray_tube_usefull_1(5*(k5-1)+5,:) = Ray_tube(5*(j-1)+5,:);
               k5 = k5 + 1;
           end
       end
end
figure
handle = patch('Faces',p,'Vertices',v,'FaceVertexCData',triangle_C1,'FaceColor','interp');
set(handle,'EdgeColor','black','LineWidth',0.1);
colorbar
view(3);
hold on
% scatter3(Ray_tube(:,1),Ray_tube(:,2),Ray_tube(:,3),3,'black','filled');
scatter3(Ray_tube_edgepoint_1(:,1),Ray_tube_edgepoint_1(:,2),Ray_tube_edgepoint_1(:,3),3,'black','filled');
scatter3(Ray_tube_usefull_1(:,1),Ray_tube_usefull_1(:,2),Ray_tube_usefull_1(:,3),3,'red','filled');
scatter3(Ray_tube_usefull_1(:,5),Ray_tube_usefull_1(:,6),Ray_tube_usefull_1(:,7),3,'black','filled');
hold off
% *********************
figure
handle = patch('Faces',p,'Vertices',v,'FaceVertexCData',triangle_C1,'FaceColor','interp');
set(handle,'EdgeColor','black','LineWidth',0.1);
colorbar
view(3);
hold on
scatter3(Ray_tube_usefull_1(:,5),Ray_tube_usefull_1(:,6),Ray_tube_usefull_1(:,7),3,'black','filled');
hold off
%% �����߹ܷ�����⡿
% ����ˮƽ����������ֱ������
E_ki = null(Ki_n);
E_ki_vv = E_ki(:,1).';% theta����
E_ki_hh = E_ki(:,2).';% phi����

% ԭʼ��ǿ�����߹����
E_field = 1;
E_area = deta_lamda^2;
% ���뵼��ķ���ϵ��
R_vv = -1;
R_hh = 1;
% ��һ����Ч���߹ܷ������
Ray_tube_usefull_1_number = length(Ray_tube_usefull_1(:,1))/5; %���߹�����
Ray_tube_usefull_1_reflect = zeros(length(Ray_tube_usefull_1(:,1)),3);%���߷��䷽��
Ray_tube_usefull_1_reflect_v = zeros(length(Ray_tube_usefull_1(:,1)),3);%���ߺ�Ŀ���ཻ��
% �����ߺ��������ߵ�·��׷��
for i = 1:length(Ray_tube_usefull_1(:,1))
    Ray_tube_usefull_1_reflect(i,:) = Ki_n - 2*(bright_triangle_n(Ray_tube_usefull_1(i,4),1)*Ki_n(1,1)+bright_triangle_n(Ray_tube_usefull_1(i,4),2)*Ki_n(1,2)+bright_triangle_n(Ray_tube_usefull_1(i,4),3)*Ki_n(1,3)).*bright_triangle_n(Ray_tube_usefull_1(i,4),:);
    Ray_tube_usefull_1_reflect_v(i,:) = Ray_tube_usefull_1(i,5:7);
end
%% ���������߳�ǿ׷�١�
% ����������ߵ��������Ԫ�ϴ�ֱ������ˮƽ�����ϵĵ糡���������������ߵ�ɢ�䳡
Ray_tube_usefull_1_reflect_e_vv = zeros(Ray_tube_usefull_1_number,3); %��ֱ����
Ray_tube_usefull_1_reflect_e_hh = zeros(Ray_tube_usefull_1_number,3); %ˮƽ����
triangle_V_E =  zeros(Ray_tube_usefull_1_number,3);
triangle_V_Es =  zeros(Ray_tube_usefull_1_number,3);
triangle_V_Hs =  zeros(Ray_tube_usefull_1_number,3);
triangle_V_R =  zeros(Ray_tube_usefull_1_number,1);
E_vv = zeros(Ray_tube_usefull_1_number,1);
E_hh = zeros(Ray_tube_usefull_1_number,1);
for j =1:Ray_tube_usefull_1_number
    Ray_tube_usefull_1_reflect_e_vv(j,:) = cross(Ki_n,bright_triangle_n(Ray_tube_usefull_1(5*(j-1)+5,4),:))/norm(cross(Ki_n,bright_triangle_n(Ray_tube_usefull_1(5*(j-1)+5,4),:)));
    Ray_tube_usefull_1_reflect_e_hh(j,:) = cross(Ki_n,Ray_tube_usefull_1_reflect_e_vv(j,:))/norm(cross(Ki_n,Ray_tube_usefull_1_reflect_e_vv(j,:)));
    triangle_V_R(j) = sqrt((Ray_tube_usefull_1(5*(j-1)+5,1) - Ray_tube_usefull_1(5*(j-1)+5,5))^2 + (Ray_tube_usefull_1(5*(j-1)+5,2) - Ray_tube_usefull_1(5*(j-1)+5,6))^2 +(Ray_tube_usefull_1(5*(j-1)+5,3) - Ray_tube_usefull_1(5*(j-1)+5,7))^2);
    triangle_V_E(j,:) = exp(-1i*Z0*triangle_V_R(j)).*E_ki_vv;
    E_vv(j) = dot(Ray_tube_usefull_1_reflect_e_vv(j,:),triangle_V_E(j,:));
    E_hh(j) = dot(Ray_tube_usefull_1_reflect_e_hh(j,:),triangle_V_E(j,:));
    triangle_V_Es(j,:) = R_vv*E_vv(j)*Ray_tube_usefull_1_reflect_e_vv(j,:) + R_hh*E_hh(j)*Ray_tube_usefull_1_reflect_e_hh(j,:);
    triangle_V_Hs(j,:) = -R_hh*E_hh(j)*Ray_tube_usefull_1_reflect_e_vv(j,:) + R_vv*E_vv(j)*Ray_tube_usefull_1_reflect_e_hh(j,:);
end
%****************************************
%% ��ͨ������εĸ���Ҷ�任����λ����Sת��Ϊ�ߵ�������͡�
% ������
% Ray_tube_usefull_1_reflect_e_vv��Ray_tube_usefull_1_reflect_e_hhΪ�������깹���ֲ�����ϵ
% �ֲ�����ϵ������������Ŀ���ཻ��Ϊ����ԭ��
% ������ԭ������ͽǵ�������Ŀ�꽻�����ӣ�����������ھֲ�����ϵ�ϵ�ͶӰ���룬
% �Ӷ�ȷ�����ǵ����߽ǵ��ھֲ�����ϵ�е�����
% Local_coordinatesΪ�����߽����ھֲ�����ϵ�е�����
Local_coordinates = zeros(Ray_tube_usefull_1_number*4,2); 
for j =1:Ray_tube_usefull_1_number
    P1(1,:) = Ray_tube_usefull_1_reflect_v(5*(j-1)+1,:) - Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:);
    P2(1,:) = Ray_tube_usefull_1_reflect_v(5*(j-1)+2,:) - Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:);
    P3(1,:) = Ray_tube_usefull_1_reflect_v(5*(j-1)+3,:) - Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:);
    P4(1,:) = Ray_tube_usefull_1_reflect_v(5*(j-1)+4,:) - Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:);
    Local_coordinates(4*(j-1)+1,1)  = dot(P1,Ray_tube_usefull_1_reflect_e_hh(j,:));%�ֲ�����ϵ������
    Local_coordinates(4*(j-1)+1,2)  = dot(P1,Ray_tube_usefull_1_reflect_e_vv(j,:));%�ֲ�����ϵ������
    Local_coordinates(4*(j-1)+2,1)  = dot(P2,Ray_tube_usefull_1_reflect_e_hh(j,:));
    Local_coordinates(4*(j-1)+2,2)  = dot(P2,Ray_tube_usefull_1_reflect_e_vv(j,:));
    Local_coordinates(4*(j-1)+3,1)  = dot(P3,Ray_tube_usefull_1_reflect_e_hh(j,:));
    Local_coordinates(4*(j-1)+3,2)  = dot(P3,Ray_tube_usefull_1_reflect_e_vv(j,:));
    Local_coordinates(4*(j-1)+4,1)  = dot(P4,Ray_tube_usefull_1_reflect_e_hh(j,:));
    Local_coordinates(4*(j-1)+4,2)  = dot(P4,Ray_tube_usefull_1_reflect_e_vv(j,:));
end
% ����S��λ���ӻ��ֺ���
% S = I14-I12-I23-I34,�����ػ���ת��Ϊ���ػ���
S_I = zeros(Ray_tube_usefull_1_number,1);
S_I14 = zeros(Ray_tube_usefull_1_number,1);
S_I12 = zeros(Ray_tube_usefull_1_number,1);
S_I23 = zeros(Ray_tube_usefull_1_number,1);
S_I34 = zeros(Ray_tube_usefull_1_number,1);
% for j =  1:Ray_tube_usefull_1_number
%    S_u(j,:) = K*dot((Ks_n -Ray_tube_usefull_1_reflect(5*(j-1)+5,:)),Ray_tube_usefull_1_reflect_e_hh(j,:));
%    S_v(j,:) = K*dot((Ks_n -Ray_tube_usefull_1_reflect(5*(j-1)+5,:)),Ray_tube_usefull_1_reflect_e_vv(j,:));
%    S_w14(j,:) = 0.5*S_u(j,:)*(Local_coordinates(4*(j-1)+1,1)-Local_coordinates(4*(j-1)+4,1))+0.5*S_v(j,:)*(Local_coordinates(4*(j-1)+1,2)-Local_coordinates(4*(j-1)+4,2));
%    if(abs(S_w14(j,:))<=0.001&&S_v(j,:)~=0)
%        T_s1 = exp(1i*(0.5*S_u(j,:)*(Local_coordinates(4*(j-1)+1,1)+Local_coordinates(4*(j-1)+4,1))+0.5*S_v(j,:)*(Local_coordinates(4*(j-1)+1,2)+Local_coordinates(4*(j-1)+4,2))));
%        S_I14(j,:) = 1i*(Local_coordinates(4*(j-1)+4,1)-Local_coordinates(4*(j-1)+1,1))/S_v(j,:)*T_s1;
%        
%    elseif(S_u(j,:)~=0&&abs(S_v(j,:))<=0.001)
%        T_k = (Local_coordinates(4*(j-1)+1,1)-Local_coordinates(4*(j-1)+4,1))/(Local_coordinates(4*(j-1)+1,2)-Local_coordinates(4*(j-1)+4,2)); 
%        T_s2 = exp(1i*S_w14(j,:))*(T_k-1i*S_u(j,:)*Local_coordinates(4*(j-1)+4,2))- exp(-1i*S_w14(j,:))*(T_k-1i*S_u(j,:)*Local_coordinates(4*(j-1)+1,2));
%        S_I14(j,:) = T_s2/S_u(j,:)^2*exp(1i*S_u(j,:)*(Local_coordinates(4*(j-1)+1,1)+Local_coordinates(4*(j-1)+4,1))*0.5);
%    elseif(S_u(j,:)==0&&S_v(j,:)==0)
%        S_I14(j,:) = (Local_coordinates(4*(j-1)+4,1)-Local_coordinates(4*(j-1)+1,1))*(Local_coordinates(4*(j-1)+1,2)-Local_coordinates(4*(j-1)+4,2))*0.5;
%    else
%        T_s1 = exp(1i*(0.5*S_u(j,:)*(Local_coordinates(4*(j-1)+1,1)+Local_coordinates(4*(j-1)+4,1))+0.5*S_v(j,:)*(Local_coordinates(4*(j-1)+1,2)+Local_coordinates(4*(j-1)+4,2))));
%        S_I14(j,:) = 1i*(Local_coordinates(4*(j-1)+4,1)-Local_coordinates(4*(j-1)+1,1))*sin(S_w14(j,:))/(S_v(j,:)*S_w14(j,:))*T_s1;
%    end
% end
for j =  1:Ray_tube_usefull_1_number
    S_I14(j,:) = S_calculating(K,Ks_n,Ray_tube_usefull_1_reflect(5*(j-1)+5,:),Ray_tube_usefull_1_reflect_e_hh(j,:),Ray_tube_usefull_1_reflect_e_vv(j,:),Local_coordinates(4*(j-1)+1,:),Local_coordinates(4*(j-1)+4,:));
    S_I12(j,:) = S_calculating(K,Ks_n,Ray_tube_usefull_1_reflect(5*(j-1)+5,:),Ray_tube_usefull_1_reflect_e_hh(j,:),Ray_tube_usefull_1_reflect_e_vv(j,:),Local_coordinates(4*(j-1)+1,:),Local_coordinates(4*(j-1)+2,:));
    S_I23(j,:) = S_calculating(K,Ks_n,Ray_tube_usefull_1_reflect(5*(j-1)+5,:),Ray_tube_usefull_1_reflect_e_hh(j,:),Ray_tube_usefull_1_reflect_e_vv(j,:),Local_coordinates(4*(j-1)+2,:),Local_coordinates(4*(j-1)+3,:));
    S_I34(j,:) = S_calculating(K,Ks_n,Ray_tube_usefull_1_reflect(5*(j-1)+5,:),Ray_tube_usefull_1_reflect_e_hh(j,:),Ray_tube_usefull_1_reflect_e_vv(j,:),Local_coordinates(4*(j-1)+3,:),Local_coordinates(4*(j-1)+4,:));
    S_I(j,:) = S_I14(j,:)-S_I12(j,:)-S_I23(j,:)-S_I34(j,:);
end
%% ����������߹ܵĳ�ǿ����A_theta��A_phi��
% E_ki_vv = E_ki(:,1).';% theta����
% E_ki_hh = E_ki(:,2).';% phi����
A_theta = zeros(Ray_tube_usefull_1_number,1);
A_phi = zeros(Ray_tube_usefull_1_number,1);
for j =  1:Ray_tube_usefull_1_number
    T_s1 = 1i*K/(4*pi)*exp(1i*K*dot(Ks_n,Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:)));
    A_theta(j,:) = dot(T_s1*(cross(-E_ki_hh,triangle_V_Es(j,:)) + Z0*cross(E_ki_vv,triangle_V_Hs(j,:))),bright_triangle_n(Ray_tube_usefull_1(5*(j-1)+5,4),:))*S_I(j,:);
    A_phi(j,:) = dot(T_s1*(cross(-E_ki_hh,triangle_V_Es(j,:)) + Z0*cross(E_ki_vv,triangle_V_Hs(j,:))),bright_triangle_n(Ray_tube_usefull_1(5*(j-1)+5,4),:))*S_I(j,:);
end
RCS_hh = 10*log10(4*pi*sum(abs(A_theta)));
RCS_vv = 4*pi*sum(abs(A_phi));
















