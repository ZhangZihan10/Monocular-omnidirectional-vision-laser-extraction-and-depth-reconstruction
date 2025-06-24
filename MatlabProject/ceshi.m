function [VBX,VBY,VRX,VRY,VGX,VGY]=ceshi(img)
load trainedNet.mat;
%[file,path]=uigetfile('D:\桌面文件夹\robot course\arduino\视觉识别\语义分割测试\测试图片\');
%filepath=fullfile(path,file);
%I=imread(filepath);
%img=imread('test1.jpg');
I=img;
figure(1);
imshow(I);

I=imresize(I,[1080, 1920]);

C=semanticseg(I,net,'MiniBatchSize', 32);
%pxds =pixelLabelDatastore(I,classes,labelIDs);

classes=["Bei","Red", "Green","Black"];
cmap=camvidColorMap;%需要更改内参数
B=labeloverlay(I,C,'ColorMap',cmap,'Transparency',0.4);
figure(2);
imshow(B),title("Semantic segmentation Result");
pixelLabelColorbar(cmap,classes);

%寻找黑色方块
C1=cellstr(C);
LB = [];
for i =7:631 %1:size(C1, 1)
    for j = 749:1250  %1:size(C1, 2)
        if strcmp(C1(i, j), "Black")%选定检测对象
            LB = [LB; i, j];
        end
    end
end

%黑色方块中心点
meanValueB = mean(LB, 1);
%图片尺寸信息为列*行
Valuex=[meanValueB(2)-70,meanValueB(2)+70];%确定区域范围
Valuey=[meanValueB(1)-70,meanValueB(1)+70];
VBX=Valuex;
VBY=Valuey;

%寻找红色方块
C1=cellstr(C);
LR = [];
for i =7:631 %1:size(C1, 1)
    for j = 749:1250  %1:size(C1, 2)
        if strcmp(C1(i, j), "Red")%选定检测对象
            LR = [LR; i, j];
        end
    end
end

%红色方块中心点
meanValueR = mean(LR, 1);
%图片尺寸信息为列*行
Valuex=[meanValueR(2)-70,meanValueR(2)+70];%确定区域范围
Valuey=[meanValueR(1)-70,meanValueR(1)+70];
VRX=Valuex;
VRY=Valuey;

%寻找绿色方块
C1=cellstr(C);
LG = [];
for i =7:631 %1:size(C1, 1)
    for j = 749:1250  %1:size(C1, 2)
        if strcmp(C1(i, j), "Green")%选定检测对象
            LG = [LG; i, j];
        end
    end
end

%绿色方块中心点
meanValueG = mean(LG, 1);
%图片尺寸信息为列*行
Valuex=[meanValueG(2)-50,meanValueG(2)+50];%确定区域范围
Valuey=[meanValueG(1)-50,meanValueG(1)+50];
VGX=Valuex;
VGY=Valuey;
end