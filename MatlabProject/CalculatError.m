%误差计算
x1=22.3698;
y1=19.8802;
x2=22.4259;
y2=19.9154;

x0=22.14;
y0=19.96;

Tnew=sqrt((x0-x1)^2+(y0-y1)^2);
Told=sqrt((x0-x2)^2+(y0-y2)^2);
disp(Tnew);
disp(Told);