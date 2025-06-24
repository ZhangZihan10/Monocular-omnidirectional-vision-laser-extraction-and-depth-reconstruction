function [centers] = steger_centerline(img, sigma, thresh)
% STEGER_CENTERLINE 基于Steger方法提取线状结构光中心
% 输入：
%   img    - 输入图像（灰度）
%   sigma  - 高斯平滑系数（典型值1.5-3.0）
%   thresh - 响应阈值（典型值0.5-1.5）
% 输出：
%   centers - 亚像素级中心线坐标(Nx2矩阵)

% 转换为双精度浮点计算
img = im2double(img);

% ========== 步骤1：高斯平滑及导数计算 ==========
% 生成高斯核及其导数
hsize = ceil(3*sigma);  % 滤波器半宽
[x,y] = meshgrid(-hsize:hsize, -hsize:hsize);

% 高斯一阶导数
Gx = -x/(2*pi*sigma^4) .* exp(-(x.^2 + y.^2)/(2*sigma^2));
Gy = -y/(2*pi*sigma^4) .* exp(-(x.^2 + y.^2)/(2*sigma^2));

% 高斯二阶导数
Gxx = (x.^2 - sigma^2)/(2*pi*sigma^6) .* exp(-(x.^2 + y.^2)/(2*sigma^2));
Gxy = (x.*y)/(2*pi*sigma^6) .* exp(-(x.^2 + y.^2)/(2*sigma^2));
Gyy = (y.^2 - sigma^2)/(2*pi*sigma^6) .* exp(-(x.^2 + y.^2)/(2*sigma^2));

% 卷积计算导数
Ix  = imfilter(img, Gx, 'conv', 'replicate');
Iy  = imfilter(img, Gy, 'conv', 'replicate');
Ixx = imfilter(img, Gxx, 'conv', 'replicate');
Ixy = imfilter(img, Gxy, 'conv', 'replicate');
Iyy = imfilter(img, Gyy, 'conv', 'replicate');

% ========== 步骤2：Hessian矩阵特征分析 ==========
[height, width] = size(img);
centers = [];

for i = 2:height-1
    for j = 2:width-1
        % 构建Hessian矩阵
        H = [Ixx(i,j) Ixy(i,j); 
             Ixy(i,j) Iyy(i,j)];
        
        % 计算特征值和特征向量
        [V, D] = eig(H);
        lambda1 = D(1,1);
        lambda2 = D(2,2);
        
        % 特征值排序 |λ1| < |λ2|
        if abs(lambda1) > abs(lambda2)
            tmp = lambda1; lambda1 = lambda2; lambda2 = tmp;
            V = fliplr(V);
        end
        
        % 线状结构条件判断
        if (lambda2 < 0) && (abs(lambda2) > thresh) && (abs(lambda1/lambda2) < 0.5)
            % ===== 步骤3：亚像素定位 =====
            nx = V(1,2);  % 法线方向
            ny = V(2,2);
            
            % 计算极值点偏移量
            t = -(Ix(i,j)*nx + Iy(i,j)*ny) / ...
                (Ixx(i,j)*nx^2 + 2*Ixy(i,j)*nx*ny + Iyy(i,j)*ny^2 + eps);
            
            % 偏移量有效性检查
            if abs(t*nx) <= 1 && abs(t*ny) <= 1
                x_hat = j + t*nx;  % 亚像素X坐标
                y_hat = i + t*ny;  % 亚像素Y坐标
                
                % 边界检查
                if x_hat > 1 && x_hat < width && y_hat > 1 && y_hat < height
                    centers = [centers; x_hat y_hat];
                end
            end
        end
    end
end

% ========== 可视化结果 ==========
figure;
imshow(img); 
hold on;
plot(centers(:,1), centers(:,2), 'r.', 'MarkerSize', 10);
title('Steger方法提取的中心线');
end