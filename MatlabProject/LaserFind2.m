%使用中心寻找法,直接使用HSV图像计算

function img=LaserFind2(image)


imwrite(image,'1.jpg');

%第一种算法，改进型,检测低光域下的激光
black_image=image;%ceshiBnew(image);%识别黑色方块，提取黑色方块所在区域.使用CNN
%imshow(black_image);

%图像处理，增强激光线段
% 将图像从RGB转换到HSV
hsvImage = rgb2hsv(black_image);
%figure;
%imshow(hsvImage);
% 红色激光的HSV范围，使用一些偏移量
hue = 0.95;
saturation = 0.95;
value = 0.95;

% 设置HSV范围（根据需要调整偏移量）
hue_offset = 0.05;
sat_offset = 0.9;
val_offset = 0.05;

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
minPixelCount = 50;   % 保留的最小尺寸
maxPixelCount = 1000; % 保留的最大尺寸

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
%%figure;
%imshow(cleanedImage);
%title('Medium light Image');%title('Medium light Image');

%原算法
%figure;
%img = las_segm(black_image);


% 第二种方法,增加一类搜索
% 红色激光的HSV范围，使用一些偏移量
hue2 = 0.95;
saturation2 = 0.95;
value2 = 0.95;

% 设置HSV范围（根据需要调整偏移量）
hue_offset2 = 0.05;
sat_offset2 = 0.05;
val_offset2 = 0.9;

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
maxPixelCount2 = 400; % 保留的最大尺寸

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
%title('Low light Image');%title('Medium light Image');


%第三种方法
%第三种方法，查找纯红色
hue3 = 0;
saturation3 = 1;
value3 = 0.6;

% 设置HSV范围（根据需要调整偏移量）
hue_offset3 = 0;
sat_offset3 = 0;
val_offset3 = 0.59;

lower_red3 = [hue3 - hue_offset3, max(saturation3 - sat_offset3, 0), max(value3 - val_offset3, 0)];
upper_red3 = [hue3 + hue_offset3, 1, 1];

% 创建掩膜，提取红色部分
redMask3 = (hsvImage(:,:,1) >= lower_red3(1)) & (hsvImage(:,:,1) <= upper_red3(1)) & ...
          (hsvImage(:,:,2) >= lower_red3(2)) & (hsvImage(:,:,2) <= upper_red3(2)) & ...
          (hsvImage(:,:,3) >= lower_red3(3)) & (hsvImage(:,:,3) <= upper_red3(3));

% 对掩膜进行形态学操作（膨胀和腐蚀）
redMask3 = imdilate(redMask3, se);
redMask3 = imerode(redMask3, se);

% 连通组件分析
cc3 = bwconncomp(redMask3);

% 计算每个连通组件的像素数量
numPixels3 = cellfun(@numel, cc3.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount3 = 40;   % 保留的最小尺寸
maxPixelCount3 = 200; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents1 = numPixels3 >= minPixelCount3 & numPixels3 <= maxPixelCount3;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage1 = false(size(redMask3));
for i = 1:length(largeComponents1)
    if largeComponents1(i)
        cleanedImage1(cc3.PixelIdxList{i}) = true;
    end
end

% 显示处理后的图像
%figure;
%imshow(cleanedImage1);
%title('Pure Red light Image');%title('Medium light Image');


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


% 使用形态学操作获得骨架化的图像
skeletonImage = bwmorph(combinedImage, 'skel', Inf);
%imshow(skeletonImage)
[y, x] = find(skeletonImage);
points = [x, y];
% 提取骨架线条的坐标
% 使用DBSCAN算法进行聚类，距离阈值设为100
epsilon = 10; % 距离阈值
minPts = 2; % 每个簇最少包含的点数
labels = dbscan(points, epsilon, minPts);

% 获取不同类别的点
uniqueLabels = unique(labels);
% 显示结果
figure;
imshow(black_image);
hold on;
% 遍历每个类别，连接每类内的点
for i = 1:length(uniqueLabels)
    if uniqueLabels(i) == -1
        continue; % 忽略噪声点
    end
    
    % 获取当前类别的点
    clusterPoints = points(labels == uniqueLabels(i), :);
    
    % 计算点之间的距离
    D = pdist2(clusterPoints, clusterPoints);
    
    % 创建图并添加边
    G = graph();
    numPoints = size(clusterPoints, 1);
    for j = 1:numPoints
        for k = j+1:numPoints
            if D(j, k) < 10
                G = addedge(G, j, k, D(j, k));
            end
        end
    end
    
    % 获取最小生成树
    T = minspantree(G);
    
    % 获取最小生成树的边
    edges = table2array(T.Edges);
    
    % 绘制骨架线条的中心线
    for j = 1:size(edges, 1)
        plot([clusterPoints(edges(j, 1), 1), clusterPoints(edges(j, 2), 1)], ...
             [clusterPoints(edges(j, 1), 2), clusterPoints(edges(j, 2), 2)], 'b-', 'LineWidth', 2);
    end
end

title('Centerline of the Laser');
hold off;


% 创建与 black_image 相同大小的二值图像
binaryImage = ones(size(black_image));
% 将RGB图像转换为灰度图像
grayImage = rgb2gray(binaryImage);
% 将灰度图像转换为二值图像
threshold = 2; % 设定阈值
binaryImage = imbinarize(grayImage, threshold);
% 在二值图像上绘制平均线段
% 遍历每个类别，连接每类内的点
for i = 1:length(uniqueLabels)
    if uniqueLabels(i) == -1
        continue; % 忽略噪声点
    end
    
    % 获取当前类别的点
    clusterPoints = points(labels == uniqueLabels(i), :);
    
    % 计算点之间的距离
    D = pdist2(clusterPoints, clusterPoints);
    
    % 创建图并添加边
    G = graph();
    numPoints = size(clusterPoints, 1);
    for j = 1:numPoints
        for k = j+1:numPoints
            if D(j, k) < 10
                G = addedge(G, j, k, D(j, k));
            end
        end
    end
    
    % 获取最小生成树
    T = minspantree(G);
    
    % 获取最小生成树的边
    edges = table2array(T.Edges);
    
    % 绘制骨架线条的中心线
    for j = 1:size(edges, 1)
        binaryImage=drawLine1(binaryImage,clusterPoints(edges(j, 1), 1), clusterPoints(edges(j, 1), 2), ...
             clusterPoints(edges(j, 2), 1), clusterPoints(edges(j, 2), 2));
    end
end

%figure;
%imshow(binaryImage);
%title('Binary Image with Average Lines');


img=binaryImage;

end