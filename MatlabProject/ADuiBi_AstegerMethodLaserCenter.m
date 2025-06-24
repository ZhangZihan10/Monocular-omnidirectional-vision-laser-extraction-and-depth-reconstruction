clear
% 参数设置
%sigma = 1;        % 匹配激光条带宽度
%threshold = 0.2;  % 根据噪声水平调整

% 读取图像

[img, cmap] = imread('002.png');

%img = ind2gray(img, cmap);
%img1=img;
%img= imread('0.jpg');
%[img, ~, alpha] = imread('005.png');imshow(img, []);
% 执行中心线提取
%centers = steger_centerline(img, sigma, threshold);


% 基本用法（全自动模式）

direction = 'vertical'; % 激光条带方向'horizontal'; %
%direction = 'horizontal'; 
sigma = 5;           % 高斯滤波系数
threshold = 0.1;       % 灰度阈值

% 读取图像

% 执行中心线提取
center_points = gray_centroid_centerline(img, direction, sigma, threshold);

% 可选：多项式拟合平滑
if size(center_points,1) > 3
    x = center_points(:,1);
    y = center_points(:,2);
    p = polyfit(x, y, 3); % 三次多项式拟合
    y_fit = polyval(p, x);
    
    figure;
    imshow(img);
    hold on;
    plot(x, y_fit, 'g-', 'LineWidth', 2);
    title('拟合后的中心线');
end

% 参数设置
sigma = 5;    % 高斯尺度（控制线宽适应）
thresh = 0.6;   % 响应阈值（控制噪声抑制）

% 执行中心线提取
centers = steger_centerline(img1, sigma, thresh);

centerLine = extractLaserCenterline_Extrema(img);

