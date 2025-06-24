function centerline = gray_centroid_centerline(image, direction, sigma, threshold)
% GRAY_CENTROID_CENTERLINE 灰度重心法提取激光中心线
% 输入参数：
%   img: 输入图像（支持彩色/灰度）
%   direction: 处理方向 ('horizontal'或'vertical')
%   sigma: 高斯滤波系数（推荐1-3）
%   threshold: 灰度阈值（0-1，相对值）
% 输出：
%   centerline: 亚像素中心线坐标[Nx2]

% 图像预处理
img = image;
if size(img,3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end
img_gray = im2double(img_gray);

% 高斯滤波去噪
if sigma > 0
    img_gray = imgaussfilt(img_gray, sigma);
end

% 自动对比度拉伸
img_gray = imadjust(img_gray);

% 初始化参数
[rows, cols] = size(img_gray);
centerline = zeros(0,2);

% 处理方向判断
is_vertical = strcmpi(direction, 'vertical');

% 主处理循环
if is_vertical
    % 垂直方向处理（逐列计算）
    for col = 1:cols
        column_data = img_gray(:, col);
        [centroid, valid] = calculate_centroid(column_data, threshold);
        if valid
            centerline(end+1, :) = [col, centroid]; %#ok<AGROW>
        end
    end
else
    % 水平方向处理（逐行计算）
    for row = 1:rows
        row_data = img_gray(row, :);
        [centroid, valid] = calculate_centroid(row_data', threshold);
        if valid
            centerline(end+1, :) = [centroid, row]; %#ok<AGROW>
        end
    end
end

% 可视化结果
figure;
imshow(img);
hold on;
plot(centerline(:,1), centerline(:,2), 'r.', 'MarkerSize', 16);
title('灰度重心法中心线提取');
hold off;

% 子函数：计算单列/行的重心
    function [centroid, valid] = calculate_centroid(data, thresh)
        max_val = max(data);
        if max_val < thresh
            centroid = 0;
            valid = false;
            return;
        end
        
        % 数据归一化
        data = data - min(data);
        data = data / max(data);
        
        % 计算加权重心
        indices = 1:length(data);
        weighted_sum = sum(indices(:) .* data(:));
        total_weight = sum(data(:));
        
        if total_weight == 0
            centroid = 0;
            valid = false;
        else
            centroid = weighted_sum / total_weight;
            valid = true;
        end
    end
end