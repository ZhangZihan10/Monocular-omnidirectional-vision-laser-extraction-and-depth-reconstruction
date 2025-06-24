% 霍夫变换，   二值图像绘制线段函数
function img = drawLine(img, point1, point2)
    % 获取线段的坐标
    x = [point1(1), point2(1)];
    y = [point1(2), point2(2)];
    
    % 使用插值法绘制线段
    numPoints = max(abs(diff(x)), abs(diff(y))) + 1;
    xIndices = round(linspace(x(1), x(2), numPoints));
    yIndices = round(linspace(y(1), y(2), numPoints));
    
    % 设置线段处的像素值为 0
    for i = 1:length(xIndices)
        img(yIndices(i), xIndices(i)) = 1;
    end
end