function centerLine = extractLaserCenterline_Extrema(imagePath)
    % 读取灰度图像
    img = imagePath;
    if size(img, 3) == 3
        img = rgb2gray(img); % 转换为灰度图像
    end
    img = double(img); % 转换为双精度数据
    
    % 计算一阶梯度
    [Gx, Gy] = gradient(img);
    gradientMagnitude = sqrt(Gx.^2 + Gy.^2);
    
    % 计算二阶梯度
    [Gxx, Gxy] = gradient(Gx);
    [Gyx, Gyy] = gradient(Gy);
    secondDerivative = Gxx + Gyy;
    
    % 查找极值点（局部最小值）
    [rows, cols] = size(img);
    centerPoints = [];
    
    for r = 2:rows-1
        for c = 2:cols-1
            if secondDerivative(r, c) < secondDerivative(r-1, c) && ...
               secondDerivative(r, c) < secondDerivative(r+1, c) && ...
               secondDerivative(r, c) < secondDerivative(r, c-1) && ...
               secondDerivative(r, c) < secondDerivative(r, c+1)
                
                centerPoints = [centerPoints; c, r];
            end
        end
    end
    
    % 显示结果
    figure; imshow(img, []); hold on;
    plot(centerPoints(:, 1), centerPoints(:, 2), 'r.', 'MarkerSize', 1);
    title('极值法激光条带中心线提取结果');
    hold off;
    
    centerLine = centerPoints;
end
