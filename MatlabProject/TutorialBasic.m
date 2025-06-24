clc
clear all
% Laser Segmentation
image = imread('TestImages/image8.jpg');
img = las_segm(image);
% Mapping
load('Omni_Calib_Results_Real.mat'); % Calib parameters
ocam_model = calib_data.ocam_model; % Calib parameters
x = 0; % Laser Plane parameters
y = 0; % Laser Plane parameters
las_dist = 200; % Laser Plane parameters
CVsyst_x = 0; % CV System position
CVsyst_y = 0; % CV System position
[x1,y1] = mapping(img,x,y,las_dist,ocam_model); % mapping function
% Finally figure:
figure;
scatter(x1,y1,5,'filled'); % Laser intersections
hold on;
plot(CVsyst_x,-CVsyst_y,'r*'); % CV System location
grid on;
% Results validation
% Left Cube
i=[600;800]; % working image region - column
j=[300;500]; % working image region - row
[C_left] = cube_dist(img,i,j,x,y,las_dist,ocam_model);
C_left = mean(C_left(:,1));
% Up Cube
i=[900;1000]; % working image region - column
j=[100;540]; % working image region - row
[C_Up] = cube_dist(img,i,j,x,y,las_dist,ocam_model);
C_Up = mean(C_Up(:,2));
% Right Cube
i=[1175;1272]; % working image region - column
j=[300;563]; % working image region - row
[C_Right] = cube_dist(img,i,j,x,y,las_dist,ocam_model);
C_Right = mean(C_Right(:,1)); 