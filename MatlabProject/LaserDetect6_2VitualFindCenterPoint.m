%不进行物体识别提取，直接进行全局图像种激光识别
clear
%cam=webcam(2);
%preview(cam);
name = "Matlab";
Client = TCPInit('127.0.0.1',55012,name);
%Client=a;
image =ImageReadTCP_One(Client,'Center'); %image=imread('image.png');

imwrite(image,'2.jpg');

%第一种算法，改进型,检测低光域下的激光
black_image=image;%ceshiBnew(image);%识别黑色方块，提取黑色方块所在区域.使用CNN
%imshow(black_image);

%图像处理，增强激光线段
% 将图像从RGB转换到HSV
hsvImage = rgb2hsv(black_image);

% 转换红色激光的RGB值到HSV值（手动估算）
red_rgb = [100/255, 50/255, 30/255]; % RGB值转换到0-1范围
red_hsv = rgb2hsv(reshape(red_rgb, [1 1 3])); % 转换到HSV

% 红色激光的HSV范围，使用一些偏移量
hue = red_hsv(1);
saturation = red_hsv(2);
value = red_hsv(3);

% 设置HSV范围（根据需要调整偏移量）
hue_offset = 0.2;
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
minPixelCount = 3;   % 保留的最小尺寸
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
figure;
imshow(cleanedImage);
title('Low light Image');

%原算法
figure;
img = las_segm(black_image);


% 第二种方法,增加一类搜索
red_rgb2 = [200/255, 5/255, 30/255]; % 新的红色RGB值
red_hsv2 = rgb2hsv(reshape(red_rgb2, [1 1 3])); % 转换到HSV

% 红色激光的HSV范围，使用一些偏移量
hue2 = red_hsv2(1);
saturation2 = red_hsv2(2);
value2 = red_hsv2(3);

% 设置HSV范围（根据需要调整偏移量）
hue_offset2 = 0.1;
sat_offset2 = 0.2;
val_offset2 = 0.2;

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
minPixelCount2 = 3;   % 保留的最小尺寸
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
figure;
imshow(cleanedImage2);
title('Medium light Image');


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
minPixelCount = 3;   % 保留的最小尺寸
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
figure;
imshow(cleanedImage1);
title('hight light Image');

% 叠加二值图像
combinedImage = cleanedImage | cleanedImage1 | cleanedImage2;
% 显示叠加后的图像
figure;
imshow(combinedImage);
title('Combined Binary Image');


% 使用形态学操作获得骨架化的图像
skeletonImage = bwmorph(combinedImage, 'skel', Inf);
% 提取骨架线条的坐标
[y, x] = find(skeletonImage);
points = [x, y];
% 使用DBSCAN算法进行聚类，距离阈值设为100
epsilon = 50; % 距离阈值
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

title('Centerline of the White Line');
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

% 显示结果
figure;
imshow(binaryImage);
title('Binary Image with Average Lines');




% Mapping
load('Omni_Calib_Results_Unity.mat'); % Calib parameters
ocam_model = calib_data.ocam_model; % Calib parameters
x = 0; % Laser Plane parameters
y = 0; % Laser Plane parameters
las_dist = 950; % Laser Plane parameters
CVsyst_x = 0; % CV System initial position 在unity中为CVSystemOrigin的位置参数z*1000
CVsyst_y = 0; % CV System initial position 在unity中为CVSystemOrigin的位置参数x*-1000
[x1,y1] = mapping(binaryImage,x,y,las_dist,ocam_model); % mapping function
% Finally figure:
figure;
scatter(x1,y1,5,'filled'); % Laser intersections
hold on;
plot(CVsyst_x,-CVsyst_y,'r*'); % CV System location
grid on;
i=[1;900];%i=[850;1250]; % working image region - column
j=[500;1900];%j=[10;700]; % working image region - row

[C_Up] = cube_dist(binaryImage,i,j,x,y,las_dist,ocam_model);
C_Up1 = mean(C_Up(:,1));%调用x的信息
C_Up2 = mean(C_Up(:,2));%调用y的信息
