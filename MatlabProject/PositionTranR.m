function  output=PositionTranR(x2,y2)
%A1=[0 0 -1 11;0 1 0 3.81; 1 0 0 3];%cv系统转unity世界坐标系，第3列为CV系统的坐标系
A2=[-1 0 0 0.175;0 -1 0 0.13; 0 0 0 1];%相机坐标系转换为世界坐标系，最后一列为相机坐标（x,y,z),单位m
%cv=[-4.49;-0.51;-7.31;1];

cv=[x2/1000;y2/1000;-0.32;1];
%cv=[-7.31;-0.32;4.48;1];
uv=A2*cv;

output=uv;%单位m
end