%寻找最小值算法,先进行深度学习识别，并提取对应物体，再进行激光识别和位置计算
%clc
%clear all
% Laser Segmentation
%image = imread('TestImages/image16.jpg');
function output=TutorialBasic2Rnew
Cub_l=10;%木块的边长一半单位mm
cam=webcam(2);
%preview(cam);

cam.Resolution='1920x1080';
%cam.Brightness=0;%调整相机亮度

image =snapshot(cam);

load('Omni_Calib_Results1.mat'); % Calib parameters
ocam_model = calib_data.ocam_model; % Calib parameters

% Results validation
% Results validation
% Left Cube
%i=[20;850]; % working image region - column
%j=[10;700]; % working image region - row
%[C_left] = cube_dist(img,i,j,x,y,las_dist,ocam_model);
%C_left = mean(C_left(:,1));
% Red Cube
black_image=ceshiRnew(image);%识别黑色方块，提取黑色方块所在区域.使用CNN
%imshow(black_image);

%图像处理，增强激光线段
% 将图像从RGB转换到HSV
hsvImage = rgb2hsv(black_image);

% 转换红色激光的RGB值到HSV值（手动估算）
red_rgb = [70/255, 0, 0]; % RGB值转换到0-1范围
red_hsv = rgb2hsv(reshape(red_rgb, [1 1 3])); % 转换到HSV

% 红色激光的HSV范围，使用一些偏移量
hue = red_hsv(1);
saturation = red_hsv(2);
value = red_hsv(3);

% 设置HSV范围（根据需要调整偏移量）
hue_offset = 0.05;
sat_offset = 0.5;
val_offset = 0.5;

lower_red1 = [hue - hue_offset, max(saturation - sat_offset, 0), max(value - val_offset, 0)];
upper_red1 = [hue + hue_offset, 1, 1];

% 创建掩膜，提取红色部分
redMask = (hsvImage(:,:,1) >= lower_red1(1)) & (hsvImage(:,:,1) <= upper_red1(1)) & ...
          (hsvImage(:,:,2) >= lower_red1(2)) & (hsvImage(:,:,2) <= upper_red1(2)) & ...
          (hsvImage(:,:,3) >= lower_red1(3)) & (hsvImage(:,:,3) <= upper_red1(3));

% 显示红色部分掩膜
figure;
imshow(redMask);
title('Red Mask');

% 对掩膜进行形态学操作（膨胀和腐蚀）
se = strel('disk', 1); % 使用较小的结构元素
redMask = imdilate(redMask, se);
redMask = imerode(redMask, se);

% 显示形态学处理后的掩膜
figure;
imshow(redMask);
title('Morphologically Processed Red Mask');

%区域激光识别
img = redMask;

%位置计算
x = 0; % Laser Plane parameters
y = 0; % Laser Plane parameters
las_dist = 131; % Laser Plane parameters
CVsyst_x = 0; % CV System position
CVsyst_y = 0; % CV System position
[x1,y1] = mapping(img,x,y,las_dist,ocam_model); % mapping function
% Finally figure:
figure(3);
scatter(x1,y1,5,'filled'); % Laser intersections
hold on;
plot(CVsyst_x,-CVsyst_y,'r*'); % CV System location
grid on;
i=[10;1250]; % working image region - column
j=[10;1000]; % working image region - row
[C_Up] = cube_dist(img,i,j,x,y,las_dist,ocam_model);
C_Up1 = min(C_Up(:,1));%调用x的信息
C_Up2 = min(C_Up(:,2));%调用y的信息
%C_Up = mean(C_Up(:,2));%调用y的信息
% Right Cube
%i=[850;1250]; % working image region - column
%j=[10;700]; % working image region - row
%[C_Right] = cube_dist(img,i,j,x,y,las_dist,ocam_model);
%C_Right1 = min(C_Right(:,1)); %调用x的信息
%C_Right = mean(C_Right(:,1)); %调用x的信息
%figure(4);
%imshow(image);
output=[C_Up1+Cub_l,C_Up2+Cub_l];%在计算值基础上y坐标加1厘米
end