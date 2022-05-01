function  S_I14 = S_calculating(K,Ks_n,Ray_tube_usefull_1_reflect,Ray_tube_usefull_1_reflect_e_hh,Ray_tube_usefull_1_reflect_e_vv,Local_coordinates_p,Local_coordinates_q)
   S_u = K*dot((Ks_n -Ray_tube_usefull_1_reflect),Ray_tube_usefull_1_reflect_e_hh);
   S_v = K*dot((Ks_n -Ray_tube_usefull_1_reflect),Ray_tube_usefull_1_reflect_e_vv);
   S_w14= 0.5*S_u*(Local_coordinates_q(1,1)-Local_coordinates_p(1,1))+0.5*S_v*(Local_coordinates_q(1,2)-Local_coordinates_p(1,2));
   if(abs(S_w14)<=0.001&&S_v~=0)
       
       T_s1 = exp(1i*(0.5*S_u*(Local_coordinates_q(1,1)+Local_coordinates_p(1,1))+0.5*S_v*(Local_coordinates_q(1,2)+Local_coordinates_p(1,2))));
       S_I14 = 1i*(Local_coordinates_q(1,1)-Local_coordinates_p(1,1))/S_v*T_s1;
       
   elseif(S_u~=0&&abs(S_v)<=0.001)
       
       T_k = (Local_coordinates_q(1,1)-Local_coordinates_p(1,1))/(Local_coordinates_q(1,2)-Local_coordinates_p(1,2)); 
       T_s2 = exp(1i*S_w14)*(T_k-1i*S_u*Local_coordinates_q(1,2))- exp(-1i*S_w14)*(T_k-1i*S_u*Local_coordinates_p(1,2));
       S_I14 = T_s2/S_u^2*exp(1i*S_u*(Local_coordinates_q(1,1)+Local_coordinates_p(1,1))*0.5);
       
   elseif(S_u==0&&S_v==0)
       
       S_I14 = (Local_coordinates_q(1,1)-Local_coordinates_p(1,1))*(Local_coordinates_q(1,2)-Local_coordinates_p(1,2))*0.5;
% 以上为三种奇异情况，以下为一般情况
   else
       
       T_s1 = exp(1i*(0.5*S_u*(Local_coordinates_q(1,1)+Local_coordinates_p(1,1))+0.5*S_v*(Local_coordinates_q(1,2)+Local_coordinates_p(1,2))));
       S_I14 = 1i*(Local_coordinates_q(1,1)-Local_coordinates_p(1,1))*sin(S_w14)/(S_v*S_w14)*T_s1;
   end

end