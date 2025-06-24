%% 自适应阈值激光图像二值化程序
% 日期：2023-10-15
% 功能：基于自适应阈值的激光条带二值化提取
clear
% 读取图像
I = imread('001.png');  % 替换为你的图像文件路径

% 若为彩色图像，转换为灰度图
if size(I, 3) == 3
    Igray = rgb2gray(I);
else
    Igray = I;
end

% 局部自适应阈值计算
% sensitivity 可调整（范围 0~1，越大越敏感，默认 0.5~0.7 比较常用）
T = adaptthresh(Igray, 0.4, 'NeighborhoodSize', 21, 'Statistic', 'median');

% 二值化图像
BW = imbinarize(Igray, T);
figure
imshow(BW);
% 可选：进行形态学操作，去除小区域或连接条纹
BW = bwareaopen(BW, 40);  % 去除小于30像素的连通区域（可根据需要调整）
% 获取连通区域信息
CC = bwconncomp(BW);
stats = regionprops(CC, 'Area');

% 找出大于1000像素的区域索引
largeIdx = find([stats.Area] > 6000);

% 构造掩码用于去除这些区域
BW_clean = BW;
for i = 1:length(largeIdx)
    BW_clean(CC.PixelIdxList{largeIdx(i)}) = 0;
end

% 显示结果
figure;
imshow(BW_clean);
title('提取的激光区域（二值图，已去除大区域）');

% 创建金黄色标记层（RGB: [1, 0.8431, 0]）
Icolor = im2double(I);  % 转为 double 用于叠加显示

overlay = Icolor;
gold = [1, 0.8431, 0];  % 金黄色

for k = 1:3
    channel = Icolor(:,:,k);
    channel(BW_clean) = gold(k);  % 仅对激光区域进行着色
    overlay(:,:,k) = channel;
end

% 显示结果
figure;
imshow(overlay);
title('原图中以金黄色标出的激光区域');