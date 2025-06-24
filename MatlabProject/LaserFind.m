%使用霍夫变换寻找直线，根据RGB进行计算

function img=LaserFind(image)


imwrite(image,'1.jpg');

%第一种算法，改进型,检测低光域下的激光
black_image=image;%ceshiBnew(image);%识别黑色方块，提取黑色方块所在区域.使用CNN
%imshow(black_image);

%图像处理，增强激光线段
% 将图像从RGB转换到HSV
hsvImage = rgb2hsv(black_image);

% 转换红色激光的RGB值到HSV值（手动估算）
red_rgb = [100/255, 5/255, 30/255]; % RGB值转换到0-1范围
red_hsv = rgb2hsv(reshape(red_rgb, [1 1 3])); % 转换到HSV

% 红色激光的HSV范围，使用一些偏移量
hue = red_hsv(1);
saturation = red_hsv(2);
value = red_hsv(3);

% 设置HSV范围（根据需要调整偏移量）
hue_offset = 0.1;
sat_offset = 0.5;
val_offset = 0.5;

lower_red1 = [hue - hue_offset, max(saturation - sat_offset, 0), max(value - val_offset, 0)];
upper_red1 = [hue + hue_offset, 1, 1];

% 创建掩膜，提取红色部分
redMask = (hsvImage(:,:,1) >= lower_red1(1)) & (hsvImage(:,:,1) <= upper_red1(1)) & ...
          (hsvImage(:,:,2) >= lower_red1(2)) & (hsvImage(:,:,2) <= upper_red1(2)) & ...
          (hsvImage(:,:,3) >= lower_red1(3)) & (hsvImage(:,:,3) <= upper_red1(3));

% 对掩膜进行形态学操作（膨胀和腐蚀）
se = strel('disk', 3); % 使用较小的结构元素
redMask = imdilate(redMask, se);
redMask = imerode(redMask, se);

% 显示形态学处理后的掩膜
% 连通组件分析
cc = bwconncomp(redMask);

% 计算每个连通组件的像素数量
numPixels = cellfun(@numel, cc.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount = 80;   % 保留的最小尺寸
maxPixelCount = 600; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents = numPixels >= minPixelCount & numPixels <= maxPixelCount;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage = false(size(redMask));
for i = 1:length(largeComponents)
    if largeComponents(i)
        cleanedImage(cc.PixelIdxList{i}) = true;
    end
end
% 显示处理后的图像
%figure;
%imshow(cleanedImage);
%title('Low light Image');

% 第二种方法,增加一类搜索
red_rgb2 = [200/255, 5/255, 30/255]; % 新的红色RGB值
red_hsv2 = rgb2hsv(reshape(red_rgb2, [1 1 3])); % 转换到HSV

% 红色激光的HSV范围，使用一些偏移量
hue2 = red_hsv2(1);
saturation2 = red_hsv2(2);
value2 = red_hsv2(3);

% 设置HSV范围（根据需要调整偏移量）
hue_offset2 = 0.1;
sat_offset2 = 0.5;
val_offset2 = 0.5;

lower_red2 = [hue2 - hue_offset2, max(saturation2 - sat_offset2, 0), max(value2 - val_offset2, 0)];
upper_red2 = [hue2 + hue_offset2, 1, 1];

% 创建掩膜，提取红色部分
redMask2 = (hsvImage(:,:,1) >= lower_red2(1)) & (hsvImage(:,:,1) <= upper_red2(1)) & ...
          (hsvImage(:,:,2) >= lower_red2(2)) & (hsvImage(:,:,2) <= upper_red2(2)) & ...
          (hsvImage(:,:,3) >= lower_red2(3)) & (hsvImage(:,:,3) <= upper_red2(3));

% 对掩膜进行形态学操作（膨胀和腐蚀）
redMask2 = imdilate(redMask2, se);
redMask2 = imerode(redMask2, se);

% 连通组件分析
cc2 = bwconncomp(redMask2);

% 计算每个连通组件的像素数量
numPixels2 = cellfun(@numel, cc2.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount2 = 50;   % 保留的最小尺寸
maxPixelCount2 = 500; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents2 = numPixels2 >= minPixelCount2 & numPixels2 <= maxPixelCount2;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage2 = false(size(redMask2));
for i = 1:length(largeComponents2)
    if largeComponents2(i)
        cleanedImage2(cc2.PixelIdxList{i}) = true;
    end
end

% 显示处理后的图像
%figure;
%imshow(cleanedImage2);
%title('Medium light Image');


%第三种方法
% 提取红色通道img=black_image;
redChannel = black_image(:,:,1);

% 将红色通道转换为灰度图像
grayRed = mat2gray(redChannel);

% 二值化处理
% 将灰度值为1的部分设置为1，其他部分为0
binaryRed = grayRed==1;
%binaryRed = imbinarize(grayRed, 'adaptive', 'Sensitivity', 0);
% 形态学操作，去除噪点和增强线条
se = strel('disk', 3); % 使用较小的结构元素
%se = strel('line', 20, 0);
dilatedRed = imdilate(binaryRed, se);
erodedRed = imerode(dilatedRed, se);
% 连通组件分析
cc = bwconncomp(erodedRed);

% 计算每个连通组件的像素数量
numPixels = cellfun(@numel, cc.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount = 15;   % 保留的最小尺寸
maxPixelCount = 500; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents = numPixels >= minPixelCount & numPixels <= maxPixelCount;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage1 = false(size(erodedRed));
for i = 1:length(largeComponents)
    if largeComponents(i)
        cleanedImage1(cc.PixelIdxList{i}) = true;
    end
end

% 显示处理后的图像
%figure;
%imshow(cleanedImage1);
%title('hight light Image');

% 叠加二值图像
combinedImage = cleanedImage | cleanedImage1 | cleanedImage2;
% 显示叠加后的图像
%figure;
%imshow(combinedImage);
%title('Combined Binary Image');

% 形态学操作，去除块状物
%se = strel('disk', 1); % 使用较大的结构元素去除块状物
%cleanedImage = imopen(combinedImage, se);

% 显示去除块状物后的图像
%figure;
%imshow(cleanedImage);
%title('Cleaned Image');

cleanedImage=combinedImage;
% 提取线状物体
% 使用霍夫变换提取线状物体
[H, theta, rho] = hough(cleanedImage);
peaks = houghpeaks(H, 200, 'threshold', ceil(0.2 * max(H(:))));
lines = houghlines(cleanedImage, theta, rho, peaks, 'FillGap', 10, 'MinLength', 40);


% 初始化起点数组
startPoints = zeros(length(lines), 2);
for k = 1:length(lines)
    startPoints(k, :) = lines(k).point1;
end

% 使用 DBSCAN 聚类
epsilon = 30; % 距离阈值
minPts = 1; % 每个簇最少包含的点数
labels = dbscan(startPoints, epsilon, minPts);

% 计算每个聚类的平均线段
uniqueLabels = unique(labels);
numClusters = length(uniqueLabels(uniqueLabels > 0));
averageLines = struct('point1', {}, 'point2', {});
for i = 1:numClusters
    clusterLines = lines(labels == uniqueLabels(i));
    if isempty(clusterLines)
        continue;
    end
    
    % 计算聚类中线段的平均起点和终点
    avgPoint1 = mean(cat(1, clusterLines.point1), 1);
    avgPoint2 = mean(cat(1, clusterLines.point2), 1);
    
    averageLines(i).point1 = avgPoint1;
    averageLines(i).point2 = avgPoint2;
end

% 显示结果
figure;
imshow(black_image);
hold on;
for k = 1:length(averageLines)
    xy = [averageLines(k).point1; averageLines(k).point2];
    plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'red');

    % 绘制每条线的起点和终点
    plot(xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
    plot(xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'green');
end
hold off;
title('Average Lines in Clusters');

% 创建与 black_image 相同大小的二值图像
binaryImage = ones(size(black_image));
% 将RGB图像转换为灰度图像
grayImage = rgb2gray(binaryImage);
% 将灰度图像转换为二值图像
threshold = 2; % 设定阈值
binaryImage = imbinarize(grayImage, threshold);
% 在二值图像上绘制平均线段
for k = 1:length(averageLines)
    % 获取线段起点和终点
    point1 = averageLines(k).point1;
    point2 = averageLines(k).point2;
    
    % 使用二值图像绘制函数绘制线段
    binaryImage = drawLine(binaryImage, point1, point2);
end

% 显示结果
figure;
imshow(binaryImage);
title('Binary Image with Average Lines');

img=binaryImage;

end