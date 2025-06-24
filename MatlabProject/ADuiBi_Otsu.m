%%全局阈值法提取激光
% 
% 清理环境
clc; clear; close all;

%% 1. 读取图像并转为灰度图
img = imread('实验1原图.png');      % 读取图像     % 使用内置图像或替换为你的图像路径
if size(img, 3) == 3
    img = rgb2gray(img);       % 转为灰度图
end
img = double(img);             % 转为双精度浮点数
[rows, cols] = size(img);      % 获取图像尺寸

%% 2. 计算灰度直方图
histogram = zeros(256, 1);    % 初始化直方图
for i = 1:rows
    for j = 1:cols
        gray_value = img(i, j) + 1;  % 灰度值从1到256（避免索引为0）
        histogram(gray_value) = histogram(gray_value) + 1;
    end
end
total_pixels = rows * cols;    % 图像总像素数

%% 3. Otsu算法：计算最佳阈值
max_variance = 0;             % 最大类间方差
optimal_threshold = 128;        % 最佳阈值

% 遍历所有可能的阈值（0~255）
for t = 1:255
    % 计算前景和背景的像素数占比
    w0 = sum(histogram(1:t)) / total_pixels;      % 背景权重
    w1 = sum(histogram(t+1:256)) / total_pixels;  % 前景权重
    
    % 计算前景和背景的灰度均值
    mu0 = sum((0:t-1)' .* histogram(1:t)) / (w0 * total_pixels);
    mu1 = sum((t:255)' .* histogram(t+1:256)) / (w1 * total_pixels);
    
    % 计算类间方差
    variance = w0 * w1 * (mu1 - mu0)^2;
    
    % 更新最佳阈值
    if variance > max_variance
        max_variance = variance;
        optimal_threshold = t - 1;  % 阈值转为0~255范围
    end
end

%% 4. 应用阈值进行二值化
binary_img = zeros(rows, cols);
binary_img(img > optimal_threshold) = 255;  % 大于阈值为前景（白色）

%% 5. 显示结果
figure;
subplot(1, 3, 1); imshow(uint8(img), []); title('原图');
subplot(1, 3, 2); imhist(uint8(img)); title('直方图'); 
subplot(1, 3, 3); imshow(uint8(binary_img), []); 
title(['Otsu二值化结果（阈值=', num2str(optimal_threshold), ')']);