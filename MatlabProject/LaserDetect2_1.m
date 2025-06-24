% 读取图像
image = imread('1.jpg');

% 转换为灰度图像
grayImage = rgb2gray(image);

% 直方图均衡化
equalizedImage = histeq(grayImage);

% 自适应直方图均衡化
adaptEqualizedImage = adapthisteq(equalizedImage);

% 高斯滤波
blurredImage = imgaussfilt(adaptEqualizedImage, 2);

% 光照归一化
normalizedImage = mat2gray(blurredImage);

% 边缘检测
edges = edge(normalizedImage, 'Canny');

% 显示结果
figure;
subplot(2, 3, 1);
imshow(grayImage);
title('原始灰度图像');

subplot(2, 3, 2);
imshow(equalizedImage);
title('直方图均衡化');

subplot(2, 3, 3);
imshow(adaptEqualizedImage);
title('自适应直方图均衡化');

subplot(2, 3, 4);
imshow(blurredImage);
title('高斯滤波');

subplot(2, 3, 5);
imshow(normalizedImage);
title('光照归一化');

subplot(2, 3, 6);
imshow(edges);
title('边缘检测');
% 转换为HSV颜色空间
hsvImage = rgb2hsv(image);

% 设置红色的颜色范围
lowerRed = [0, 0.6, 0.6];
upperRed = [0.05, 1, 1];

% 颜色过滤
redMask = (hsvImage(:,:,1) >= lowerRed(1) & hsvImage(:,:,1) <= upperRed(1)) & ...
          (hsvImage(:,:,2) >= lowerRed(2) & hsvImage(:,:,2) <= upperRed(2)) & ...
          (hsvImage(:,:,3) >= lowerRed(3) & hsvImage(:,:,3) <= upperRed(3));

% 显示颜色过滤结果
figure;
imshow(redMask);
title('红色激光线条的颜色过滤结果');
