%根据相机的旋转角度，进行目标物体坐标系的旋转，使其与机器人坐标系一致，单位毫秒
%thed为相机坐标系绕Z轴顺时针旋转到机器人坐标的角度，x1，y1为激光点在相机坐标系的位置，xc，yc为相机在机器人坐标系内位置
function [x3,y3]=trans(thed,x1,y1,xc,yc)
%x1=[-1766];y1=[-3344.4];thed=180;
x3=x1;y3=y1;
for i=1:length(x1)
%alpha=atand(y1(i)/x1(i));
%x3(i)=sqrt(x1(i)^2+y1(i)^2)*cosd(thed-alpha);
%y3(i)=sqrt(x1(i)^2+y1(i)^2)*sind(alpha-thed);

x3(i)=cosd(thed)*x1(i)+sind(thed)*y1(i)+xc;
y3(i)=-sind(thed)*x1(i)+cosd(thed)*y1(i)+yc;
end

end