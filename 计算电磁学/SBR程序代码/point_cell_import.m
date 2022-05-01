% 判断各角点相交的面元
function point_cell = point_cell_import(v,bright_triangle,point_xyz,Ki_n)
bright_triangle_number = length(bright_triangle);
parfor number_triangle = 1:bright_triangle_number
       [normal_t(number_triangle),normal_u(number_triangle),normal_v(number_triangle)] = normal_judge(v,bright_triangle(number_triangle,:),point_xyz,Ki_n); 
end
for number_triangle1 = 1:bright_triangle_number
    if(normal_u(number_triangle1)>=0&&normal_v(number_triangle1)>=0&&(normal_v(number_triangle1)+normal_u(number_triangle1))<=1)
        Reflection_point = (1-normal_u(number_triangle1)-normal_v(number_triangle1))*v(bright_triangle(number_triangle1,1),:)+normal_u(number_triangle1)*v(bright_triangle(number_triangle1,2),:)+normal_v(number_triangle1)*v(bright_triangle(number_triangle1,3),:);
        point_cell = [point_xyz number_triangle1 Reflection_point];
        break;
    else
        point_cell = [point_xyz 0 0 0 0];
    end
end
end