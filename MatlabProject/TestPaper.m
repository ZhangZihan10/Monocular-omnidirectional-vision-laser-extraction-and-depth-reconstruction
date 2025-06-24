%不进行物体识别提取，直接进行全局图像种激光识别
%clear
cam=webcam(2);
%preview(cam);

cam.Resolution='1920x1080';
%cam.Brightness=20;    %调整相机亮度

image =snapshot(cam);
imwrite(image,'1.jpg');

%image=imread('002.png');
black_image=image;%ceshiBnew(image);%识别黑色方块，提取黑色方块所在区域.使用CNN
%imshow(image);

%原算法
figure;
img = las_segm(black_image);

figure;
[rows, cols] = find(img == 1); % 找到img中值为1的元素的行和列坐标
black_image1=black_image;
% 遍历所有找到的坐标
for k = 1:length(rows)
    % 将black_image中对应位置的像素值设置为黑色
    black_image1(rows(k), cols(k), :) = [50, 100, 50];
end
% 显示修改后的black_image
imshow(black_image1);

%图像处理，增强激光线段
% 将图像从RGB转换到HSV
hsvImage = rgb2hsv(black_image);
figure;
imshow(hsvImage);

%% 第一种方法 红色激光的HSV范围，使用一些偏移量
hue = 0.9;
saturation = 0.01;
value = 0.9;

% 设置HSV范围（根据需要调整偏移量）
hue_offset = 0.04;
sat_offset = 0.95;
val_offset = 0.02;

lower_red1 = [hue, saturation, value];
upper_red1 = [1, 1, 1];

% 创建掩膜，提取红色部分
redMask = (hsvImage(:,:,1) >= lower_red1(1)) & (hsvImage(:,:,1) <= upper_red1(1)) & ...
          (hsvImage(:,:,2) >= lower_red1(2)) & (hsvImage(:,:,2) <= upper_red1(2)) & ...
          (hsvImage(:,:,3) >= lower_red1(3)) & (hsvImage(:,:,3) <= upper_red1(3));

% 对掩膜进行形态学操作（膨胀和腐蚀）
se = strel('disk', 7); % 使用较小的结构元素
redMask = imdilate(redMask, se);
redMask = imerode(redMask, se);
%imshow(redMask);
% 显示形态学处理后的掩膜
% 连通组件分析
cc = bwconncomp(redMask);

% 计算每个连通组件的像素数量
numPixels = cellfun(@numel, cc.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount = 50;   % 保留的最小尺寸
maxPixelCount = 900; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents = numPixels >= minPixelCount & numPixels <= maxPixelCount;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage = false(size(redMask));
for i = 1:length(largeComponents)
    if largeComponents(i)
        cleanedImage(cc.PixelIdxList{i}) = true;
    end
end
% 显示处理后的图像
figure;
imshow(cleanedImage);
title('Red interference light Image');%title('Medium light Image');



%% 第二种方法 
% 红色激光的HSV范围，使用一些偏移量
hue2 = 0.95;
saturation2 = 0.95;
value2 = 0.95;

% 设置HSV范围（根据需要调整偏移量）
hue_offset2 = 0.05;
sat_offset2 = 0.05;
val_offset2 = 0.95;

lower_red2 = [hue2 - hue_offset2, max(saturation2 - sat_offset2, 0), max(value2 - val_offset2, 0)];
upper_red2 = [1, 1, 1];

% 创建掩膜，提取红色部分
redMask2 = (hsvImage(:,:,1) >= lower_red2(1)) & (hsvImage(:,:,1) <= upper_red2(1)) & ...
          (hsvImage(:,:,2) >= lower_red2(2)) & (hsvImage(:,:,2) <= upper_red2(2)) & ...
          (hsvImage(:,:,3) >= lower_red2(3)) & (hsvImage(:,:,3) <= upper_red2(3));

% 对掩膜进行形态学操作（膨胀和腐蚀）
redMask2 = imdilate(redMask2, se);
redMask2 = imerode(redMask2, se);

% 连通组件分析
cc2 = bwconncomp(redMask2);

% 计算每个连通组件的像素数量
numPixels2 = cellfun(@numel, cc2.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount2 = 100;   % 保留的最小尺寸
maxPixelCount2 = 300; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents2 = numPixels2 >= minPixelCount2 & numPixels2 <= maxPixelCount2;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage2 = false(size(redMask2));
for i = 1:length(largeComponents2)
    if largeComponents2(i)
        cleanedImage2(cc2.PixelIdxList{i}) = true;
    end
end

% 显示处理后的图像
figure;
imshow(cleanedImage2);
title('Low light Image');%title('Medium light Image');


%% 第三种方法，查找纯红色
hue3 = 0;
saturation3 = 0.9;
value3 = 0.1;

% 设置HSV范围（根据需要调整偏移量）
hue_offset3 = 0.01;
sat_offset3 = 0.01;
val_offset3 = 0.55;

lower_red3 = [0, max(saturation3, 0), value3];
upper_red3 = [hue_offset3, 1, 0.4];

% 创建掩膜，提取红色部分
redMask3 = (hsvImage(:,:,1) >= lower_red3(1)) & (hsvImage(:,:,1) <= upper_red3(1)) & ...
          (hsvImage(:,:,2) >= lower_red3(2)) & (hsvImage(:,:,2) <= upper_red3(2)) & ...
          (hsvImage(:,:,3) >= lower_red3(3)) & (hsvImage(:,:,3) <= upper_red3(3));

figure;
[rows, cols] = find(redMask3 == 1); % 找到img中值为1的元素的行和列坐标
black_image3=black_image;
% 遍历所有找到的坐标
for k = 1:length(rows)
    % 将black_image中对应位置的像素值设置为黑色
    black_image3(rows(k), cols(k), :) = [0, 0, 255];
end
% 显示修改后的black_image
imshow(black_image3);
figure;imshow(redMask3);
% 对掩膜进行形态学操作（膨胀和腐蚀）
redMask3 = imdilate(redMask3, se);
redMask3 = imerode(redMask3, se);
%figure;imshow(redMask3);
% 连通组件分析
cc3 = bwconncomp(redMask3);

% 计算每个连通组件的像素数量
numPixels3 = cellfun(@numel, cc3.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）imshow(redMask3)
minPixelCount3 = 100;   % 保留的最小尺寸
maxPixelCount3 = 700; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents1 = numPixels3 >= minPixelCount3 & numPixels3 <= maxPixelCount3;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage1 = false(size(redMask3));
for i = 1:length(largeComponents1)
    if largeComponents1(i)
        cleanedImage1(cc3.PixelIdxList{i}) = true;
    end
end

% 显示处理后的图像
figure;
imshow(cleanedImage1);
title('Pure Red light Image');%title('Medium light Image');