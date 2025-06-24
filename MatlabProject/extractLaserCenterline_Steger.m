function centerLine = extractLaserCenterline_Steger(imagePath)
    % 读取灰度图像
    img = imagePath;
    if size(img, 3) == 3
        img = rgb2gray(img); % 转换为灰度图像
    end
    img = double(img); % 转换为双精度数据
    
    % 计算图像梯度
    [Gx, Gy] = gradient(img);
    [Gxx, Gxy] = gradient(Gx);
    [Gyx, Gyy] = gradient(Gy);
    
    % 计算 Hessian 矩阵的特征值
    rows = size(img, 1);
    cols = size(img, 2);
    centerLine = [];
    
    for r = 1:rows
        for c = 1:cols
            H = [Gxx(r, c), Gxy(r, c); Gyx(r, c), Gyy(r, c)];
            [V, D] = eig(H);
            [~, minIdx] = min(diag(D)); % 选择最小特征值对应的特征向量
            dirVector = V(:, minIdx); 
            
            % 计算亚像素中心位置
            x_new = c + dirVector(1);
            y_new = r + dirVector(2);
            centerLine = [centerLine; x_new, y_new];
        end
    end
    
    % 显示结果
    figure; imshow(img, []); hold on;
    plot(centerLine(:, 1), centerLine(:, 2), 'r.', 'MarkerSize', 1);
    title('Steger方法激光条带中心线提取结果');
    hold off;
end
