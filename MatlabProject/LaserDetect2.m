%不进行物体识别提取，直接进行全局图像种激光识别
clear
cam=webcam(2);
%preview(cam);

cam.Resolution='1920x1080';
%cam.Brightness=-29;%调整相机亮度

image =snapshot(cam);
imwrite(image,'1.jpg');
load('Omni_Calib_Results1.mat'); % Calib parameters
ocam_model = calib_data.ocam_model; % Calib parameters

%第一种算法，改进型
black_image=image;%ceshiBnew(image);%识别黑色方块，提取黑色方块所在区域.使用CNN
%imshow(black_image);

%图像处理，增强激光线段
% 将图像从RGB转换到HSV
hsvImage = rgb2hsv(black_image);

% 转换红色激光的RGB值到HSV值（手动估算）
red_rgb = [100/255, 0, 0]; % RGB值转换到0-1范围
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
se = strel('disk', 3); % 使用较小的结构元素
redMask = imdilate(redMask, se);
redMask = imerode(redMask, se);

% 显示形态学处理后的掩膜
figure;
imshow(redMask);
title('Morphologically Processed Red Mask');
% 连通组件分析
cc = bwconncomp(redMask);

% 计算每个连通组件的像素数量
numPixels = cellfun(@numel, cc.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount = 50;   % 保留的最小尺寸
maxPixelCount = 800; % 保留的最大尺寸

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
title('Cleaned Binary Image');

%第二种算法
figure;
img = las_segm(black_image);

%第三种方法
% 提取红色通道img=black_image;
redChannel = black_image(:,:,1);

% 将红色通道转换为灰度图像
grayRed = mat2gray(redChannel);

% 二值化处理
% 将灰度值为1的部分设置为1，其他部分为0
binaryRed = grayRed==1;
%binaryRed = imbinarize(grayRed, 'adaptive', 'Sensitivity', 0);
% 形态学操作，去除噪点和增强线条
se = strel('line', 20, 0);
dilatedRed = imdilate(binaryRed, se);
erodedRed = imerode(dilatedRed, se);
figure;
imshow(erodedRed);
title('Binary Image with Specific Gray Value');
% 连通组件分析
cc = bwconncomp(erodedRed);

% 计算每个连通组件的像素数量
numPixels = cellfun(@numel, cc.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount = 15;   % 保留的最小尺寸
maxPixelCount = 500; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents = numPixels >= minPixelCount & numPixels <= maxPixelCount;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage1 = false(size(erodedRed));
for i = 1:length(largeComponents)
    if largeComponents(i)
        cleanedImage1(cc.PixelIdxList{i}) = true;
    end
end

% 显示处理后的图像
figure;
imshow(cleanedImage1);
title('Cleaned Binary Image');

% 使用Canny边缘检测
edges = edge(cleanedImage1, 'canny');

% 霍夫变换检测直线
[H, theta, rho] = hough(edges);
peaks = houghpeaks(H, 5, 'threshold', ceil(0.3*max(H(:))));
lines = houghlines(edges, theta, rho, peaks, 'FillGap', 5, 'MinLength', 10);

% 显示结果
figure, imshow(black_image), hold on
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'green');
end
hold off