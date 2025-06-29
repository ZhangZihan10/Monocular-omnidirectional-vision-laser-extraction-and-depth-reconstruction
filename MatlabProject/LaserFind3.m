%使用中心寻找法,直接使用HSV图像计算

function img=LaserFind3(image)

imwrite(image,'1.jpg');

%第一种算法，改进型,检测低光域下的激光
black_image=image;%ceshiBnew(image);%识别黑色方块，提取黑色方块所在区域.使用CNN
%imshow(black_image);

%图像处理，增强激光线段
% 将图像从RGB转换到HSV
hsvImage = rgb2hsv(black_image);
%figure;
%imshow(hsvImage);

%% 第一种方法 红色激光的HSV范围，使用一些偏移量
hue = 0.97;
saturation = 0.95;
value = 0.97;

% 设置HSV范围（根据需要调整偏移量）
hue_offset = 0.05;
sat_offset = 0.95;
val_offset = 0.05;

lower_red1 = [hue - hue_offset, max(saturation - sat_offset, 0), max(value - val_offset, 0)];
upper_red1 = [1, 1, 1];

% 创建掩膜，提取红色部分
redMask = (hsvImage(:,:,1) >= lower_red1(1)) & (hsvImage(:,:,1) <= upper_red1(1)) & ...
          (hsvImage(:,:,2) >= lower_red1(2)) & (hsvImage(:,:,2) <= upper_red1(2)) & ...
          (hsvImage(:,:,3) >= lower_red1(3)) & (hsvImage(:,:,3) <= upper_red1(3));

% 对掩膜进行形态学操作（膨胀和腐蚀）
se = strel('disk', 7); % 使用较小的结构元素
redMask = imdilate(redMask, se);
redMask = imerode(redMask, se);

% 显示形态学处理后的掩膜
% 连通组件分析
cc = bwconncomp(redMask);

% 计算每个连通组件的像素数量
numPixels = cellfun(@numel, cc.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount = 100;   % 保留的最小尺寸
maxPixelCount = 900; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents = numPixels >= minPixelCount & numPixels <= maxPixelCount;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage = false(size(redMask));
for i = 1:length(largeComponents)
    if largeComponents(i)
        cleanedImage(cc.PixelIdxList{i}) = true;
    end
end
% 显示处理后的图


%% 第二种方法 
% 红色激光的HSV范围，使用一些偏移量
hue2 = 0.95;
saturation2 = 0.95;
value2 = 0.95;

% 设置HSV范围（根据需要调整偏移量）
hue_offset2 = 0.05;
sat_offset2 = 0.05;
val_offset2 = 0.95;

lower_red2 = [hue2 - hue_offset2, max(saturation2 - sat_offset2, 0), max(value2 - val_offset2, 0)];
upper_red2 = [1, 1, 1];

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
minPixelCount2 = 60;   % 保留的最小尺寸
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



%% 第三种方法，查找纯红色
hue3 = 0;
saturation3 = 1;
value3 = 0.6;

% 设置HSV范围（根据需要调整偏移量）
hue_offset3 = 0.001;
sat_offset3 = 0.001;
val_offset3 = 0.55;

lower_red3 = [0, max(saturation3 - sat_offset3, 0), max(value3 - val_offset3, 0)];
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
minPixelCount3 = 70;   % 保留的最小尺寸
maxPixelCount3 = 500; % 保留的最大尺寸

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



%% 第四种方法
hue4 = 0.92;
saturation4 = 0.4;
value4 = 0.4;

% 设置HSV范围（根据需要调整偏移量）
hue_offset4 = 0;
sat_offset4 = 0;
val_offset4 = 0;

lower_red4 = [hue4 - hue_offset4, max(saturation4 - sat_offset4, 0), max(value4 - val_offset4, 0)];
upper_red4 = [1, 1, 1];

% 创建掩膜，提取红色部分
redMask4 = (hsvImage(:,:,1) >= lower_red4(1)) & (hsvImage(:,:,1) <= upper_red4(1)) & ...
          (hsvImage(:,:,2) >= lower_red4(2)) & (hsvImage(:,:,2) <= upper_red4(2)) & ...
          (hsvImage(:,:,3) >= lower_red4(3)) & (hsvImage(:,:,3) <= upper_red4(3));

% 对掩膜进行形态学操作（膨胀和腐蚀）
redMask4 = imdilate(redMask4, se);
redMask4 = imerode(redMask4, se);

% 连通组件分析
cc4 = bwconncomp(redMask4);

% 计算每个连通组件的像素数量
numPixels4 = cellfun(@numel, cc4.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount4 = 100;   % 保留的最小尺寸
maxPixelCount4 = 300; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents1 = numPixels4 >= minPixelCount4 & numPixels4 <= maxPixelCount4;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage4 = false(size(redMask4));
for i = 1:length(largeComponents1)
    if largeComponents1(i)
        cleanedImage4(cc4.PixelIdxList{i}) = true;
    end
end

% 显示处理后的图像



%% 第五种方法
hue5 = 0;
saturation5 = 0;
value5 = 0.9;


lower_red5 = [hue5, saturation5, value5];
upper_red5 = [0.2, 0.5, 1];

% 创建掩膜，提取红色部分
redMask5 = (hsvImage(:,:,1) >= lower_red5(1)) & (hsvImage(:,:,1) <= upper_red5(1)) & ...
          (hsvImage(:,:,2) >= lower_red5(2)) & (hsvImage(:,:,2) <= upper_red5(2)) & ...
          (hsvImage(:,:,3) >= lower_red5(3)) & (hsvImage(:,:,3) <= upper_red5(3));

% 对掩膜进行形态学操作（膨胀和腐蚀）
redMask5 = imdilate(redMask5, se);
redMask5 = imerode(redMask5, se);

% 连通组件分析
cc5 = bwconncomp(redMask5);

% 计算每个连通组件的像素数量
numPixels5 = cellfun(@numel, cc5.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount5 = 100;   % 保留的最小尺寸
maxPixelCount5 = 300; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents1 = numPixels5 >= minPixelCount5 & numPixels5 <= maxPixelCount5;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage5 = false(size(redMask5));
for i = 1:length(largeComponents1)
    if largeComponents1(i)
        cleanedImage5(cc5.PixelIdxList{i}) = true;
    end
end

% 显示处理后的图像
%% 第六种方法灰度阈值选取亮度高点
% 如果是彩色图像，转换为灰度图像
if size(image, 3) == 3
    imgg = rgb2gray(image);
end
% 搜寻图像中灰度值大于240的位置
threshold = 248;
binary_img5 = imgg > threshold;  % 大于240的位置变为1，其他为0

% 对掩膜进行形态学操作（膨胀和腐蚀）
binary_img5 = imdilate(binary_img5, se);
binary_img5 = imerode(binary_img5, se);

% 连通组件分析
cc6 = bwconncomp(binary_img5);

% 计算每个连通组件的像素数量
numPixels6 = cellfun(@numel, cc6.PixelIdxList);

% 设置要保留的最小和最大像素数量（根据需要调整）
minPixelCount6 = 80;   % 保留的最小尺寸
maxPixelCount6 = 600; % 保留的最大尺寸

% 过滤掉不符合尺寸要求的组件
largeComponents1 = numPixels6 >= minPixelCount6 & numPixels6 <= maxPixelCount6;

% 创建一个新的二值图像，只保留符合尺寸要求的组件
cleanedImage6 = false(size(binary_img5));
for i = 1:length(largeComponents1)
    if largeComponents1(i)
        cleanedImage6(cc6.PixelIdxList{i}) = true;
    end
end
%显示二值图像
%figure;
%imshow(cleanedImage6);
%title('二值图像 (灰度值大于240的点)');
% 显示灰度图像的二维曲面图
%%figure;
%surf(double(img));
%colormap gray;  % 使用灰度颜色
%shading interp;  % 插值处理，使图像更平滑
%title('2D Surface of Grayscale Image');
%xlabel('X');
%ylabel('Y');
%zlabel('Gray Value');

%% 图像组合与分析

% 叠加二值图像
combinedImage = cleanedImage | cleanedImage1 | cleanedImage2 | cleanedImage4 | cleanedImage5|cleanedImage6 ;
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
epsilon = 20; % 距离阈值
minPts = 20; % 每个簇最少包含的点数
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
            if D(j, k) < 30
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
            if D(j, k) < 30
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