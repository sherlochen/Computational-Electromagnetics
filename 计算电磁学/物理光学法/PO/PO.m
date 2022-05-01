
filename = 'C:\Users\LXW\Desktop\物理光学法\算例\Five_spheres\5qiu_wangge.stl';
F = linspace(5e9,10e9,11);
RCS = zeros(1,length(F));
h = waitbar(0,'PO计算中...');
for i = 1:length(F)
    RCS(i) =  PO_function(F(i),filename);
    waitbar(i/length(F))
end
 close(h);