%% 激光条纹中心提取改进版
% 作者：根据优化方案实现
% 日期：2023-10-15
% 功能：鲁棒的亚像素级中心线提取

%% 参数设置
clc; clear; close all;
imagePath = '005.png';       % 输入图像路径
minPeakHeight = 0.2;         % 峰高阈值（相对最大强度）
smoothWindow = 15;           % 后处理平滑窗口

%% 图像预处理
try
    % 读取并转换图像
    [img, cmap] = imread(imagePath);
    if ~isempty(cmap)
        grayImg = ind2gray(img, cmap);
    else
        grayImg = rgb2gray(img);
    end
    
    % 自适应对比度增强
    enhancedImg = adapthisteq(grayImg, 'ClipLimit',0.02);
    
catch ME
    error('图像处理失败: %s', ME.message);
end

%% 初始化参数
[rows, cols] = size(enhancedImg);
centers = zeros(1, cols);
validFlags = false(1, cols);  % 有效检测标志

%% 主处理循环（动态窗口）
for col = 1:cols
    % 获取当前列数据
    columnData = double(enhancedImg(:, col));
    
    % 动态确定窗口
    [windowStart, windowEnd] = findDynamicWindow(columnData, minPeakHeight);
    
    if windowEnd - windowStart < 5
        continue;  % 跳过无效列
    end
    
    % 提取窗口数据
    y = (windowStart:windowEnd)';
    I = columnData(y) + eps;  % 避免零值
    
    % 加权最小二乘拟合
    [y0, isValid] = robustGaussianFit(y, I);
    
    if isValid
        centers(col) = y0;
        validFlags(col) = true;
    end
end

%% 后处理优化
% 中值滤波去噪
centers(validFlags) = medfilt1(centers(validFlags), smoothWindow);

% 线性插值补全缺失值
xValid = find(validFlags);
centers = interp1(xValid, centers(xValid), 1:cols, 'linear', 'extrap');

%% 可视化结果
figure('Position', [100 100 1200 600])
subplot(1,2,1), imshow(grayImg), title('原始图像')
subplot(1,2,2), imshow(enhancedImg)
hold on
plot(1:cols, centers, 'r-', 'LineWidth', 1.2)
title('改进中心线提取')
hold off

%% 核心函数定义
function [startIdx, endIdx] = findDynamicWindow(columnData, minPeakHeight)
    % 动态窗口检测函数
    grad = gradient(columnData);
    [maxInt, peakPos] = max(columnData);
    
    % 有效区域检测
    threshold = minPeakHeight * maxInt;
    validRegion = columnData > threshold;
    
    if ~any(validRegion)
        startIdx = 1;
        endIdx = length(columnData);
        return;
    end
    
    % 确定边界
    startIdx = find(validRegion, 1, 'first');
    endIdx = find(validRegion, 1, 'last');
    
    % 边界扩展
    margin = ceil(0.2*(endIdx - startIdx));
    startIdx = max(1, startIdx - margin);
    endIdx = min(length(columnData), endIdx + margin);
end

function [y0, isValid] = robustGaussianFit(y, I)
    % 鲁棒高斯拟合函数
    maxIter = 10;
    tol = 1e-6;
    
    % 初始化参数
    weights = sqrt(I);  % 光子噪声权重
    X = [ones(size(y)), y, y.^2];
    
    % 迭代加权最小二乘
    for iter = 1:maxIter
        % 构建加权矩阵
        W = diag(weights);
        
        % 求解系数
        beta = (X' * W * X) \ (X' * W * log(I));
        
        % 更新权重
        residual = log(I) - X*beta;
        newWeights = weights .* exp(-0.5*(residual).^2);
        
        % 检查收敛
        if norm(newWeights - weights) < tol
            break;
        end
        weights = newWeights;
    end
    
    % 计算中心点
    a1 = beta(2);
    a2 = beta(3);
    y0 = -a1/(2*a2);
    
    % 有效性验证
    isValid = (y0 >= min(y)) && (y0 <= max(y)) && (abs(a2) > 1e-6);
end