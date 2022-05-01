function  [normal_t,normal_u,normal_v] = normal_judge(v,p,point_xyz,Ki_n)
e1 = v(p(1,2),:) - v(p(1,1),:);
e2 = v(p(1,3),:) - v(p(1,1),:);
s =  point_xyz(1,:) - v(p(1,1),:);
normal_t = det(s,e1,e2)/det(-Ki_n,e1,e2);
normal_u = det(-Ki_n,s,e2)/det(-Ki_n,e1,e2);
normal_v = det(-Ki_n,e1,s)/det(-Ki_n,e1,e2);
    function result = det(a,b,c)
           result = dot(a,cross(b,c));
    end
end