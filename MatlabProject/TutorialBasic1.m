%平均算法
%clc
%clear all
% Laser Segmentation
%image = imread('TestImages/image16.jpg');
function output=TutorialBasic1

cam=webcam(2);
%preview(cam);

cam.Resolution='1920x1080';

image =snapshot(cam);


img = las_segm(image);
% Mapping
load('Omni_Calib_Results1.mat'); % Calib parameters
ocam_model = calib_data.ocam_model; % Calib parameters
x = 0; % Laser Plane parameters
y = 0; % Laser Plane parameters
las_dist = 145; % Laser Plane parameters
CVsyst_x = 0; % CV System position
CVsyst_y = 0; % CV System position
[x1,y1] = mapping(img,x,y,las_dist,ocam_model); % mapping function
% Finally figure:
figure(3);
scatter(x1,y1,5,'filled'); % Laser intersections
hold on;
plot(CVsyst_x,-CVsyst_y,'r*'); % CV System location
grid on;
% Results validation
% Results validation
% Left Cube
%i=[500;900]; % working image region - column
%j=[100;700]; % working image region - row
%[C_left] = cube_dist(img,i,j,x,y,las_dist,ocam_model);
%C_left = mean(C_left(:,1));
% Up Cube
%i=[850;1250]; % working image region - column
%j=[10;700]; % working image region - row
[VX,VY]=ceshi(image);%识别黑色方块，划定黑色方块所在区域.使用CNN
i=[VX(1);VX(2)]; % working image region - column
j=[VY(1);VY(2)]; % working image region - row
%将数组变为整数
i = round(i);
j = round(j);

[C_Up] = cube_dist(img,i,j,x,y,las_dist,ocam_model);
C_Up = mean(C_Up(:,2));%调用y的信息
% Right Cube
i=[850;1250]; % working image region - column
j=[10;700]; % working image region - row
[C_Right] = cube_dist(img,i,j,x,y,las_dist,ocam_model);
C_Right = mean(C_Right(:,1)); %调用x的信息
figure(4);
imshow(image);
output=[C_Right,C_Up+10];%在计算值基础上y坐标加1厘米
end