clear all
% close all
filename = "F:\计算电磁学\SBR程序代码\5qiu_wangge.stl";
%% 【网格】预处理
[v, p, n] = stlread(filename);
% 入射波方向设置
theta = 90/180*pi;
phi = 0/180*pi;
Ki_x = sin(theta)*cos(phi);
Ki_y = sin(theta)*sin(phi);
Ki_z = cos(theta);
Ki_n = [Ki_x,Ki_y,Ki_z]/norm([Ki_x,Ki_y,Ki_z]);
Ks_n = -Ki_n;
% *************************************
triangle_number = length(p); %面元数量
% 判断面元是否被遮挡
triangle_O1 = zeros(triangle_number,1);  %面元判断遮挡符号 1为亮区、0为暗区
triangle_C1 = zeros(3*triangle_number,1);
k1 = 1;
for i = 1:triangle_number
    % ********************************************************************
    % 单遮挡处理 将入射波方向同面元法向量相乘 判断亮区和暗区
    triangle_c = n(i,1)*Ki_n(1,1) + n(i,2)*Ki_n(1,2) + n(i,3)*Ki_n(1,3);
    if(triangle_c>0)     %入射波同三角面元法向量成锐角，该面元赋值为0被遮挡
        triangle_O1(i,1)= 0;
        triangle_C1(1+(i-1)*3:3+(i-1)*3,1)=0;
    elseif(triangle_c<0) %入射波同三角面元法向量成钝角，该面元赋值为1未遮挡
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
%% 【投影】
%确定射线入射口面
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
%% 【划分入射口面网格】
c0 = 3e8;
f = 10e9;%linspace(1e9,2e9,11);
lamda = c0./f;
deta_lamda = lamda/10;
Z0 = 120*pi; %自由空间波阻抗
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
%% 【射线和面元相交判断】
point_cell{M,N} = zeros;
for i = 1:M
    for j = 1:N
%         % 将入射网格中各角点保存
          point_xyz = (inv(F)*[Distance_r(i,j);Y_theta(i);Z_phi(j)]).' ;
%         % 判断各角点相交的面元
          point_cell{i,j} = point_cell_import(v,bright_triangle,point_xyz,Ki_n);
    end
end
% ****************************************
%% 【射线管归类】
% 五个坐标构成一组射线管，前四个坐标为逆时针变化的射线管角点，最后一个坐标为中心点
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
%% 【判断射线管的有效性】
% Ray_tube_usefull_1 和 Ray_tube 每行有7位数组成，前三个为入射点位置，第四位为相交面元，后三位为反射点位置
Ray_tube_number = length(Ray_tube(:,1))/5; %射线管数量
k5 = 1; %有效射线管数量
for j =1:Ray_tube_number
       w1 = Ray_tube(5*(j-1)+1,4); % 角点射线对应着的面元
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
               % 由于射线管追踪深度不一致，造成的无线射线管暂不考虑
               Ray_tube_usefull_1(5*(k5-1)+1,:) = Ray_tube(5*(j-1)+1,:); %有效射线管
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
%% 【射线管反射求解】
% 定义水平极化波、垂直极化波
E_ki = null(Ki_n);
E_ki_vv = E_ki(:,1).';% theta方向
E_ki_hh = E_ki(:,2).';% phi方向

% 原始场强及射线管面积
E_field = 1;
E_area = deta_lamda^2;
% 理想导体的反射系数
R_vv = -1;
R_hh = 1;
% 第一次有效射线管反射求解
Ray_tube_usefull_1_number = length(Ray_tube_usefull_1(:,1))/5; %射线管数量
Ray_tube_usefull_1_reflect = zeros(length(Ray_tube_usefull_1(:,1)),3);%射线反射方向
Ray_tube_usefull_1_reflect_v = zeros(length(Ray_tube_usefull_1(:,1)),3);%射线和目标相交点
% 角射线和中心射线的路径追踪
for i = 1:length(Ray_tube_usefull_1(:,1))
    Ray_tube_usefull_1_reflect(i,:) = Ki_n - 2*(bright_triangle_n(Ray_tube_usefull_1(i,4),1)*Ki_n(1,1)+bright_triangle_n(Ray_tube_usefull_1(i,4),2)*Ki_n(1,2)+bright_triangle_n(Ray_tube_usefull_1(i,4),3)*Ki_n(1,3)).*bright_triangle_n(Ray_tube_usefull_1(i,4),:);
    Ray_tube_usefull_1_reflect_v(i,:) = Ray_tube_usefull_1(i,5:7);
end
%% 【中心射线场强追踪】
% 求解中心射线到达各个面元上垂直分量和水平分量上的电场，并求解出中心射线的散射场
Ray_tube_usefull_1_reflect_e_vv = zeros(Ray_tube_usefull_1_number,3); %垂直分量
Ray_tube_usefull_1_reflect_e_hh = zeros(Ray_tube_usefull_1_number,3); %水平分量
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
%% 【通过多边形的傅里叶变换将相位因子S转换为边的线性求和】
% 首先以
% Ray_tube_usefull_1_reflect_e_vv和Ray_tube_usefull_1_reflect_e_hh为横纵坐标构建局部坐标系
% 局部坐标系以中心射线与目标相交点为中心原点
% 将中心原点坐标和角点射线与目标交点链接，计算各向量在局部坐标系上的投影距离，
% 从而确定各角点射线角点在局部坐标系中的坐标
% Local_coordinates为角射线交点在局部坐标系中的坐标
Local_coordinates = zeros(Ray_tube_usefull_1_number*4,2); 
for j =1:Ray_tube_usefull_1_number
    P1(1,:) = Ray_tube_usefull_1_reflect_v(5*(j-1)+1,:) - Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:);
    P2(1,:) = Ray_tube_usefull_1_reflect_v(5*(j-1)+2,:) - Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:);
    P3(1,:) = Ray_tube_usefull_1_reflect_v(5*(j-1)+3,:) - Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:);
    P4(1,:) = Ray_tube_usefull_1_reflect_v(5*(j-1)+4,:) - Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:);
    Local_coordinates(4*(j-1)+1,1)  = dot(P1,Ray_tube_usefull_1_reflect_e_hh(j,:));%局部坐标系横坐标
    Local_coordinates(4*(j-1)+1,2)  = dot(P1,Ray_tube_usefull_1_reflect_e_vv(j,:));%局部坐标系纵坐标
    Local_coordinates(4*(j-1)+2,1)  = dot(P2,Ray_tube_usefull_1_reflect_e_hh(j,:));
    Local_coordinates(4*(j-1)+2,2)  = dot(P2,Ray_tube_usefull_1_reflect_e_vv(j,:));
    Local_coordinates(4*(j-1)+3,1)  = dot(P3,Ray_tube_usefull_1_reflect_e_hh(j,:));
    Local_coordinates(4*(j-1)+3,2)  = dot(P3,Ray_tube_usefull_1_reflect_e_vv(j,:));
    Local_coordinates(4*(j-1)+4,1)  = dot(P4,Ray_tube_usefull_1_reflect_e_hh(j,:));
    Local_coordinates(4*(j-1)+4,2)  = dot(P4,Ray_tube_usefull_1_reflect_e_vv(j,:));
end
% 计算S相位因子积分函数
% S = I14-I12-I23-I34,将三重积分转换为二重积分
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
%% 【计算各射线管的场强函数A_theta和A_phi】
% E_ki_vv = E_ki(:,1).';% theta方向
% E_ki_hh = E_ki(:,2).';% phi方向
A_theta = zeros(Ray_tube_usefull_1_number,1);
A_phi = zeros(Ray_tube_usefull_1_number,1);
for j =  1:Ray_tube_usefull_1_number
    T_s1 = 1i*K/(4*pi)*exp(1i*K*dot(Ks_n,Ray_tube_usefull_1_reflect_v(5*(j-1)+5,:)));
    A_theta(j,:) = dot(T_s1*(cross(-E_ki_hh,triangle_V_Es(j,:)) + Z0*cross(E_ki_vv,triangle_V_Hs(j,:))),bright_triangle_n(Ray_tube_usefull_1(5*(j-1)+5,4),:))*S_I(j,:);
    A_phi(j,:) = dot(T_s1*(cross(-E_ki_hh,triangle_V_Es(j,:)) + Z0*cross(E_ki_vv,triangle_V_Hs(j,:))),bright_triangle_n(Ray_tube_usefull_1(5*(j-1)+5,4),:))*S_I(j,:);
end
RCS_hh = 10*log10(4*pi*sum(abs(A_theta)));
RCS_vv = 4*pi*sum(abs(A_phi));
















