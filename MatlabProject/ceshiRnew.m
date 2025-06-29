function [segmanticImg]=ceshiRnew(img)
load maskRCNNModel.mat;
%[file,path]=uigetfile('D:\桌面文件夹\robot course\arduino\视觉识别\语义分割测试\测试图片\');
%filepath=fullfile(path,file);
%I=imread(filepath);
%img=imread('test1.jpg');
I=img;
figure;
imshow(I);

I=imresize(I,[1080, 1920]);

C=semanticseg(I,net,'MiniBatchSize', 32);
%pxds =pixelLabelDatastore(I,classes,labelIDs);

classes=["Bei","Red", "Green","Black"];
cmap=camvidColorMap;%需要更改内参数
B=labeloverlay(I,C,'ColorMap',cmap,'Transparency',0.4);
figure;
imshow(B),title("Semantic segmentation Result");
pixelLabelColorbar(cmap,classes);

% 提取黑色方块区域掩膜
    blackMask = C == 'Red';
    
    % 初始化全黑图像
    black_image = zeros(size(I), 'uint8');
    
    % 将黑色方块区域复制到全黑图像中
    for channel = 1:size(I, 3)
        black_image(:,:,channel) = I(:,:,channel) .* uint8(blackMask);
    end
    segmanticImg=black_image;
    % 显示处理后的图像
    figure;
    imshow(black_image), title('Black Box Segmentation Result');

%寻找黑色方块
%C1=cellstr(C);
%LB = [];
%for i =7:631 %1:size(C1, 1)
%    for j = 749:1250  %1:size(C1, 2)
%        if strcmp(C1(i, j), "Red")%选定检测对象
%            LB = [LB; i, j];
%        end
%    end
%end

%黑色方块中心点
%meanValueB = mean(LB, 1);
%图片尺寸信息为列*行
%Valuex=[meanValueB(2)-70,meanValueB(2)+70];%确定区域范围
%Valuey=[meanValueB(1)-70,meanValueB(1)+70];
%VBX=Valuex;
%VBY=Valuey;
end