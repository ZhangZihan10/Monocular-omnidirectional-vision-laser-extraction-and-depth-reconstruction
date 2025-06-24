%寻找最小值算法，深度学习找到黑色物体的中心点坐标，对黑色物体进行矩形提取，之后进行激光识别处理
%clc
%clear all
% Laser Segmentation
%image = imread('TestImages/image16.jpg');
function output=TutorialBasic2BnewTest
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
[VX,VY]=ceshiB(image);%识别黑色方块，划定黑色方块所在区域.使用CNN
i=[VX(1);VX(2)]; % working image region - column
j=[VY(1);VY(2)]; % working image region - row
%将数组变为整数
i = round(i);
j = round(j);

%物体图像提取
% 获取图像尺寸
[rows, cols, channels] = size(image);
% 初始化一个全黑图像
black_image = zeros(rows, cols, channels, 'uint8');
% 处理边界情况，确保范围在图像尺寸内
i(1) = max(1, i(1));
i(2) = min(cols, i(2));
j(1) = max(1, j(1));
j(2) = min(rows, j(2));
% 将范围内的图像复制到黑色图像中
black_image(j(1):j(2), i(1):i(2), :) = image(j(1):j(2), i(1):i(2), :);
% 保存处理后的图像
imwrite(black_image, 'output_image.jpg');
% 显示处理后的图像
imshow(black_image);

%图像处理
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
sat_offset = 0.4;
val_offset = 0.4;

lower_red1 = [hue - hue_offset, max(saturation - sat_offset, 0), max(value - val_offset, 0)];
upper_red1 = [hue + hue_offset, 1, 1];

% 创建掩膜，提取红色部分
redMask = (hsvImage(:,:,1) >= lower_red1(1)) & (hsvImage(:,:,1) <= upper_red1(1)) & ...
          (hsvImage(:,:,2) >= lower_red1(2)) & (hsvImage(:,:,2) <= upper_red1(2)) & ...
          (hsvImage(:,:,3) >= lower_red1(3)) & (hsvImage(:,:,3) <= upper_red1(3));

% 对掩膜进行形态学操作（膨胀和腐蚀）
se = strel('disk', 2);
redMask = imdilate(redMask, se);
redMask = imerode(redMask, se);

% 提取红色部分
redObjectsMask = bsxfun(@times, black_image, cast(redMask, 'like', black_image));

% 将红色部分叠加到原图像
outputImage = black_image;
outputImage(repmat(redMask, [1, 1, 3])) = redObjectsMask(repmat(redMask, [1, 1, 3]));

% 显示处理后的图像
imshow(outputImage);

%区域激光识别
img = las_segm(outputImage);

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