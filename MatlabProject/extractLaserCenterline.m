function centerLine = extractLaserCenterline(imagePath)
    % 读取灰度图像
    img = imagePath;
    if size(img, 3) == 3
        img = rgb2gray(img); % 转换为灰度图像
    end
    img = double(img); % 转换为双精度数据
    
    [rows, cols] = size(img);
    centerLine = zeros(rows, 1); % 预分配中心线数组
    
    for r = 1:rows
        intensitySum = sum(img(r, :)); % 计算当前行的灰度值总和
        if intensitySum > 0
            centerLine(r) = sum((1:cols) .* img(r, :)) / intensitySum; % 计算重心
        else
            centerLine(r) = NaN; % 无激光条带时设为 NaN
        end
    end
    
    % 显示结果
    figure;imshow(img); hold on;% imshow(img, []); hold on;
    plot(centerLine, 1:rows, 'r-', 'LineWidth', 1.5); % 叠加绘制中心线
    title('激光条带中心线提取结果');
    hold off;
end