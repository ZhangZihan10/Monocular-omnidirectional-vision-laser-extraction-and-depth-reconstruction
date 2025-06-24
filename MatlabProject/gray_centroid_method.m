clear

[img, cmap] = imread('0.jpg');
img = ind2gray(img, cmap);


extractLaserCenterline(img);

%extractLaserCenterline_Steger(img);
% 读取图像
