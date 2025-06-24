% 图片路径列表（你可以添加更多路径）
image_files = {
    %'a1.jpg';'a5.jpg';'a2.jpg';'a4.jpg';'a6.jpg';
    % 'b1.jpg';'b5.jpg';'b2.jpg';'b4.jpg';'b3.jpg';
    'c1.jpg';'c5.jpg';'c2.jpg';'c4.jpg';'c3.jpg';
    % 添加更多图片路径
};

% 初始化 HSV 和颜色数据集合
H_all = [];
S_all = [];
V_all = [];
RGB_all = [];

% 遍历每张图片
for i = 1:length(image_files)
    img = imread(image_files{i});
    
    % 可选：缩小图像
    img = imresize(img, 0.25);  % 可根据需要调整
    
    % 转换为 HSV
    hsv_img = rgb2hsv(img);
    H = hsv_img(:,:,1);
    S = hsv_img(:,:,2);
    V = hsv_img(:,:,3);
    
    % 展平为向量
    H = H(:);
    S = S(:);
    V = V(:);
    
    % 展平 RGB，并归一化
    rgb_flat = reshape(img, [], 3);
    rgb_norm = double(rgb_flat) / 255;
    
    % 可选：移除背景或暗区像素
    valid_idx = V > 0.05;
    H = H(valid_idx);
    S = S(valid_idx);
    V = V(valid_idx);
    rgb_norm = rgb_norm(valid_idx, :);
    
    % 添加到总集合
    H_all = [H_all; H];
    S_all = [S_all; S];
    V_all = [V_all; V];
    RGB_all = [RGB_all; rgb_norm];
end

% 绘图
figure;
scatter3(H_all, S_all, V_all, 8, RGB_all, 'filled');
xlabel('Hue');
ylabel('Saturation');
zlabel('Value');
title('Combined HSV Distribution from Multiple Images');
xlim([0 1]); ylim([0 1]); zlim([0 1]);
grid on;
view(45, 25);
