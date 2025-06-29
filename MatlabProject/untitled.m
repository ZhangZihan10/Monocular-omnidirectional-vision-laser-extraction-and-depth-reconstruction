x = [12.3 15.52 15.75 9.11 11.94] ; 
y = [0  -8 -9.5 7 0.1];
p = polyfit(x,y,1);
disp(p);
%x2 = 0:0.001:1;
%y2 = polyval(p,x2);
plot(x,y,'o');%,'o',x2,y2)
%grid on
%s = sprintf('y = (%.1f) x^3 + (%.1f) x^2 + (%.1f) x + (%.1f)',p(1),p(2));%,p(3),p(4));
%text(2,400,s)

%0.03*19230.8-1007.7

% 创建符号变量
syms x y;

% 定义方程组
eq1 = 0.12*x + y == 13000;
eq2 = 0.057*x + y ==50;

% 求解方程组
sol = solve([eq1, eq2], [x, y]);

% 显示解
disp(sol);

% 定义方程组
eq1 =  y == 137;
eq2 = 1.57*x + y ==34;

% 求解方程组
sol = solve([eq1, eq2], [x, y]);

% 显示解
disp(sol);