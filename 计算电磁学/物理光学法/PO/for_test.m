clear
c0 = 3e8;%光速
f = 10e9;%入射波频率
lamda = c0/f;
K0 = 2*pi/(c0/f);%波数
filename = "F:\物理光学法\算例\complex_model\complex_model.stl";
[v, p, n] = stlread(filename);
figure
patch('Faces',p,'Vertices',v,'FaceColor','white')
%patch(p,v,'FaceColor','white')
view(3);
theta  =0:10;%入射角度设置theta phi
phi = 10;
for j = 1:length(theta)
    Ki_x(j) = sin(theta(j)/180*pi)*cos(phi/180*pi);
    Ki_y(j) = sin(theta(j)/180*pi)*sin(phi/180*pi);
    Ki_z(j) = cos(theta(j)/180*pi);
    Ki_n(j,1) = [Ki_x(j)]/norm([Ki_x(j),Ki_y(j),Ki_z(j)]);%入射x方向的单位矢量
    Ki_n(j,2) = [Ki_y(j)]/norm([Ki_x(j),Ki_y(j),Ki_z(j)]);%入射y方向的单位矢量
    Ki_n(j,3) = [Ki_z(j)]/norm([Ki_x(j),Ki_y(j),Ki_z(j)]);%入射z方向的单位矢量
    Ks_n = -Ki_n;%散射方向矢量
    % ********************************

    for k=1:size(Ki_n,1)
        Ki_EH{k,1}=null(Ki_n(k,:));
        Ki_E1{k,1}=(Ki_EH{k,1}(:,1))';
        Ki_E = cell2mat(Ki_E1)      %入射电场单位矢量，在PO算法中不需要考虑计划特性
    end
    % 判断面元是否被遮挡
    triangle_number = length(p); %面元数量
    triangle_O1 = zeros(triangle_number,1);  %面元判断遮挡符号 1为亮区、0为暗区
    triangle_C1 = zeros(3*triangle_number,1);
    %***************************************
    % 计算各面元的中心坐标
    triangle_Z_points = zeros(triangle_number,3);
    triangle_R = zeros(triangle_number,1);
    triangle_Rs = zeros(triangle_number,1);
    triangle_S = zeros(triangle_number,1);
    %**************************************
    % 求解感应电流
    % 计算入射波达到各面元的电场、磁场
    % A_E = 1;%入射电场幅度
    Z0 = 120*pi; %自由空间中的波阻抗
    triangle_E = zeros(triangle_number,3);  %到达各面元的电场
    triangle_H = zeros(triangle_number,3);  %到达各面元的磁场
    %**************************************
    % 在亮区感应电流 Jn=2n*Hin, 在暗区感应电流为0
    triangle_J = zeros(triangle_number,3);
    triangle_Es = zeros(triangle_number,3);
    triangle_Rcs = zeros(triangle_number,1);
    h = waitbar(0,'PO计算中...');
    for i = 1:triangle_number
        % ********************************************************************
        % 方法一：单遮挡处理 将入射波方向同面元法向量相乘
        triangle_c = n(i,1)*Ki_n(j,1) + n(i,2)*Ki_n(j,2) + n(i,3)*Ki_n(j,3);
        if(triangle_c>0)     %入射波同三角面元法向量成锐角，该面元赋值为0被遮挡
            triangle_O1(i)= 0;
            triangle_C1(1+(i-1)*3:3+(i-1)*3,:)=0;
        elseif(triangle_c<0) %入射波同三角面元法向量成钝角，该面元赋值为1未遮挡
            triangle_O1(i)= 1;
            triangle_C1(1+(i-1)*3:3+(i-1)*3,:)=1;
        end
        % ********************************************************************
        triangle_L1 = sqrt((v(p(i,2),1)-v(p(i,3),1))^2+(v(p(i,2),2)-v(p(i,3),2))^2+(v(p(i,2),3)-v(p(i,3),3))^2);
        triangle_L2 = sqrt((v(p(i,1),1)-v(p(i,3),1))^2+(v(p(i,1),2)-v(p(i,3),2))^2+(v(p(i,1),3)-v(p(i,3),3))^2);
        triangle_L3 = sqrt((v(p(i,1),1)-v(p(i,2),1))^2+(v(p(i,1),2)-v(p(i,2),2))^2+(v(p(i,1),3)-v(p(i,2),3))^2);
        triangle_Z_points(i,1)=(triangle_L1*v(p(i,1),1)+triangle_L2*v(p(i,2),1)+triangle_L3*v(p(i,3),1))/(triangle_L1+triangle_L2+triangle_L3);
        triangle_Z_points(i,2)=(triangle_L1*v(p(i,1),2)+triangle_L2*v(p(i,2),2)+triangle_L3*v(p(i,3),2))/(triangle_L1+triangle_L2+triangle_L3);
        triangle_Z_points(i,3)=(triangle_L1*v(p(i,1),3)+triangle_L2*v(p(i,2),3)+triangle_L3*v(p(i,3),3))/(triangle_L1+triangle_L2+triangle_L3);
        triangle_R(i) =  dot(Ki_n(j,:),triangle_Z_points(i,:));   % 计算各面元中心到坐标原点的距离在入射方向上的投影；r'*ki
        triangle_Rs(i) = dot(triangle_Z_points(i,:),Ks_n(j,:));  % 计算各面元中心到坐标原点的距离在散射方向上的投影
        triangle_P = (triangle_L1+triangle_L2+triangle_L3)/2;
        triangle_S(i) = sqrt(triangle_P*(triangle_P-triangle_L1)*(triangle_P-triangle_L2)*(triangle_P-triangle_L3));
        % ********************************************************************
        triangle_E(i,:) = exp(-1i*K0*triangle_R(i)).*Ki_E(j,:); %面元入射电场
        triangle_H(i,:) = cross(Ki_n(j,:),triangle_E(i,:));     %面元入射磁场
        if(triangle_O1(i)==1)
            triangle_J(i,:) = cross(n(i,:),triangle_H(i,:));%面元表面电流
            triangle_Es(i,:) = cross(triangle_J(i,:),Ks_n(j,:))*exp(1i*K0*triangle_Rs(i))*triangle_S(i); %面元远场
            %triangle_Rcs(i) = dot(n(i,:),Ks_n(j,:))*exp(2*1i*K0*dot(Ks_n(j,:),triangle_Z_points(i,:)))*triangle_S(i);
        end
        waitbar(i/triangle_number)
    end
    close(h);
    %% ******************************
    figure
    handle = patch('Faces',p,'Vertices',v,'FaceVertexCData',triangle_C1,'FaceColor','interp');
    set(handle,'EdgeColor','black','LineWidth',0.1);
    colorbar
    view(3);
    hold on
    scatter3(triangle_Z_points(:,1),triangle_Z_points(:,2),triangle_Z_points(:,3),3,'black','filled');
    hold off
    %% ******************************
    R = 10e100;
    Es_1 = sum(triangle_Es);
    ES1(j,:)= -1i*K0*exp(-1i*K0*R)/(2*pi*R)*cross(Ks_n(j,:),Es_1);
    RCS1(j,:) =4*pi*R^2*(abs(ES1(j,1))+abs(ES1(j,2))+abs(ES1(j,3)))^2;
    %% ***************************************
    %ES2 = -1i*K0*exp(-1i*K0*R)/(2*pi*R)*sum(triangle_Rcs).*Ki_E;
    %RCS2 = 20*log10(4*pi/(lamda^2)*(abs(sum(triangle_Rcs)))^2);
end
plot(theta,10*log10(RCS1),'linewidth',2);
title('f = 1GHz');
set(gca,'fontname','Times New Roman');   % 设置坐标轴字体
ylabel('\fontname{Times New Roman}Backscattering Coefficient (dB)');   %设置y轴标签
xlabel('\fontname{Times New Roman}theta(degree)');   %设置x轴标签
grid on;
