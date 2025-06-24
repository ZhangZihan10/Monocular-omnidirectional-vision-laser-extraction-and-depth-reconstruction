%% 清理环境  自适应局部阈值法
clc; clear; close all;

%% 1. 读取图像并转为灰度图
% 使用内置图像测试（例如：'text.png'），也可替换为其他图像路径
img = imread('002.png');      % 读取图像
if size(img, 3) == 3
    img = rgb2gray(img);       
end
img = double(img);             
[rows, cols] = size(img);      

%% 2. 定义参数
window_size = 91;   % 改为奇数
k = 0.1;
R = 9;

%% 3. 计算积分图（扩展零行零列避免索引越界）
% 创建扩展后的积分图（行和列各+1）
integral_img = zeros(rows+1, cols+1);
integral_sq_img = zeros(rows+1, cols+1);

% 计算累积和（从第2行第2列开始填充）
integral_img(2:end, 2:end) = cumsum(cumsum(img, 1), 2);
integral_sq_img(2:end, 2:end) = cumsum(cumsum(img.^2, 1), 2);

%% 4. 遍历像素
binary_img = zeros(rows, cols); 
half_win = floor(window_size / 2); 

for i = 1:rows
    for j = 1:cols
        % 窗口边界（基于原始图像坐标）
        y1 = max(i - half_win, 1);
        y2 = min(i + half_win, rows);
        x1 = max(j - half_win, 1);
        x2 = min(j + half_win, cols);
        
        % 转换为积分图坐标（+1偏移）
        iy1 = y1 + 1; iy2 = y2 + 1;
        ix1 = x1 + 1; ix2 = x2 + 1;
        
        % 计算sum_val（自动处理边界）
        sum_val = integral_img(iy2, ix2);
        if y1 > 1
            sum_val = sum_val - integral_img(iy1-1, ix2);
        end
        if x1 > 1
            sum_val = sum_val - integral_img(iy2, ix1-1);
        end
        if y1 > 1 && x1 > 1
            sum_val = sum_val + integral_img(iy1-1, ix1-1);
        end
        
        % 计算sum_sq（同理）
        sum_sq = integral_sq_img(iy2, ix2);
        if y1 > 1
            sum_sq = sum_sq - integral_sq_img(iy1-1, ix2);
        end
        if x1 > 1
            sum_sq = sum_sq - integral_sq_img(iy2, ix1-1);
        end
        if y1 > 1 && x1 > 1
            sum_sq = sum_sq + integral_sq_img(iy1-1, ix1-1);
        end
        
        % 计算均值和标准差
        area = (y2 - y1 + 1) * (x2 - x1 + 1);
        mean_val = sum_val / area;
        var_val = (sum_sq - sum_val^2 / area) / area;
        std_dev = sqrt(max(var_val, 0)); % 避免负方差
        
        % 计算阈值
        threshold = mean_val * (1 + k * (std_dev / R - 1));
        
        % 二值化
        if img(i, j) > threshold
            binary_img(i, j) = 255;
        end
    end
end

%% 5. 显示结果
figure;
subplot(1,2,1); imshow(uint8(img), []); title('原图');
subplot(1,2,2); imshow(uint8(binary_img), []); title('修复后结果');