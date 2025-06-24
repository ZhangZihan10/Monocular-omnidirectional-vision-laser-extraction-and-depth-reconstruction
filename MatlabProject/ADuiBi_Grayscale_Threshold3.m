%%固定灰度阈值选取亮度高点


% 读取图像
image = imread('组合1_2.jpg');  % 替换为你的图像文件路径
se = strel('disk', 7); % 使用较小的结构元素

% 如果是彩色图像，转换为灰度图像
if size(image, 3) == 3
    imgg = rgb2gray(image);
end
% 搜寻图像中灰度值大于240的位置
threshold = 240;
binary_img5 = imgg > threshold;  % 大于240的位置变为1，其他为0

% 对掩膜进行形态学操作（膨胀和腐蚀）
binary_img5 = imdilate(binary_img5, se);
binary_img5 = imerode(binary_img5, se);

% 连通组件分析
cc6 = bwconncomp(binary_img5);

% 计算每个连通组件的像素数量
numPixels6 = cellfun(@numel, cc6.PixelIdxList);

figure;
imshow(binary_img5);
% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount6 = 10;   % 保留的最小尺寸
maxPixelCount6 = 5000; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents1 = numPixels6 >= minPixelCount6 & numPixels6 <= maxPixelCount6;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage6 = false(size(binary_img5));
for i = 1:length(largeComponents1)
    if largeComponents1(i)
        cleanedImage6(cc6.PixelIdxList{i}) = true;
    end
end
%显示二值图像
figure;
imshow(cleanedImage6);
title('二值图像 (灰度值大于240的点)');