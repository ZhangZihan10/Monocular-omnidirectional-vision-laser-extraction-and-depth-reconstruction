%% 自适应阈值激光图像提取程序
% 日期：2025-05-23
% 功能：基于Sauvola算法的激光图像自适应提取


%% 激光条纹中心提取完整流程
clear; clc; close all;
% 如果是彩色图像，转换为灰度图像
img=imread("0.jpg");
% 参数设置
n = 15;               % 每侧采样点数(总点数=2n+1)
laser_threshold = 30; % 激光条纹存在判定阈值

%% 1. 图像读取与预处理
% 读取测试图像（需要替换为实际图像路径）
% 转换为灰度图像（使用R通道）
gray_img = img(:,:,1); 

% 显示原始图像
figure(1)
imshow(gray_img)
title('原始图像')

%% 2. ROI选择（可选，按实际需要取消注释）
% roi = drawrectangle('Color','r');  % 手动选择ROI
% gray_img = imcrop(gray_img, roi.Position);

%% 3. 图像增强
enhanced_img = imadjust(gray_img, [0.3 0.7], []); % 对比度扩展
figure(2)
imshow(enhanced_img)
title('增强后的图像')

%% 4. 逐列处理
[rows, cols] = size(enhanced_img);
centers = zeros(1, cols);    % 存储中心坐标
valid_cols = false(1, cols); % 有效列标记

% 进度条初始化
h = waitbar(0,'正在处理列...');

for col = 1:cols
    % 获取当前列数据
    column_data = double(enhanced_img(:, col));
    
    % 跳过低对比度列
    if max(column_data) < laser_threshold
        continue
    end
    
    % 找到最大光强位置
    [~, max_idx] = max(column_data);
    
    % 确定采样范围
    start_idx = max(1, max_idx - n);
    end_idx = min(rows, max_idx + n);
    
    % 提取采样点
    y = (start_idx:end_idx)';
    I = column_data(start_idx:end_idx);
    
    % 添加微小量避免log(0)
    I(I == 0) = 1e-6;
    
    % 高斯拟合
    H = log(I);
    Y = [ones(length(y),1), y, y.^2]; % 设计矩阵
    
    % 最小二乘求解
    coeffs = (Y'*Y) \ (Y'*H);
    
    % 计算中心位置
    a1 = coeffs(2);
    a2 = coeffs(3);
    y0 = -a1/(2*a2);
    
    % 有效性检查
    if y0 >= 1 && y0 <= rows
        centers(col) = y0;
        valid_cols(col) = true;
    end
    
    % 更新进度条
    waitbar(col/cols, h);
end

close(h) % 关闭进度条

%% 5. 结果后处理
% 线性插值填补无效列
centers(~valid_cols) = interp1(find(valid_cols), centers(valid_cols), find(~valid_cols), 'linear');

%% 6. 可视化结果
figure(3)
imshow(enhanced_img)
hold on
plot(1:cols, centers, 'r-', 'LineWidth', 1.5)
title('提取的激光中心线')
hold off

%% 7. 保存结果（可选）
% save('laser_centers.mat', 'centers');