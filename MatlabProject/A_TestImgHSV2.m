% 读取你上传的 HSV 可视化图像（注意：其实是 RGB 格式显示的）
img = imread('b1.jpg');

% 转换为 HSV 色彩空间（范围都是 [0, 1]）
hsv_img = rgb2hsv(img);

% 展平各通道
H = hsv_img(:,:,1);
S = hsv_img(:,:,2);
V = hsv_img(:,:,3);
H = H(:);
S = S(:);
V = V(:);

% 获取原图颜色作为 scatter 点颜色
rgb_flat = double(reshape(img, [], 3)) / 255;

% 可选：去除全黑背景点（Value 极低）
valid_idx = V > 0.05;  % 可以调整阈值
H = H(valid_idx);
S = S(valid_idx);
V = V(valid_idx);
rgb_flat = rgb_flat(valid_idx, :);

% 绘制 HSV 分布图
figure;
scatter3(H, S, V, 8, rgb_flat, 'filled');
xlabel('Hue');
ylabel('Saturation');
zlabel('Value');
title('HSV Distribution of Uploaded Image');
xlim([0 1]); ylim([0 1]); zlim([0 1]);
grid on;
view(45, 25);
