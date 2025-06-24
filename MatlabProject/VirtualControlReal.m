%clc;

%close all;
%name = "Matlab";
%Client = TCPInit('127.0.0.1',55016,name);

%arduino=serialport("COM13",115200);                                         %只需要运行1次，连接端口,现在用
%传输图像unity
%fig= capture_fig(Client);

%grabSend(arduino,"a");%手臂到达初始位置
%pause(3);

gripping_point = 0.056;

%gripping_point = 0.1978;

L(1) = Link([0 0 0.067 0 1]);  %定义连杆的D-H参数，关节角，连杆偏距，连杆长度，连杆转角
L(2) = Link([0 -0.017 0.092 0 0]);
L(3) = Link([0 -0.01 0.095 0 0]);
L(4) = Link([0 -0.04 0 0 0]);
L(5) = Link([0 0 0 0 0]);
L(6) = Link([0 0 0 0 0]);

L(1).qlim = [0.01 0.12];              %关节角度限制
L(2).qlim = [-105 115]/180*pi;
L(3).qlim = [-75 205]/180*pi;
L(4).qlim = [0 180]/180*pi;

robot = SerialLink(L);              %连接连杆
joints = [0.12 0 0 0 0 0];            %指定的关节角

%坐标点输入
X1 = 0.1981;  % Z
Y1 = 0.00808; % -X
Z1 = -0.01;  % Y - GREEN
t = 0:0.1:6;

T = transl(X1, Y1, Z1)* trotz(180);%根据给定终点，得到终点位姿
qi1 = robot.ikine(T,'mask',[1 1 1 1 0 0]);%根据终点点位姿，得到终点关节角
qf1 = [0.12 0 0 0 0 0];%机器人初始位置

qi1(1,4)=qi1(1,2)+qi1(1,3);%+Loc(3)*pi/180;%添加抓取点末端位姿，加入了姿态调整
q = jtraj(qf1,qi1,t);%五次多项式轨迹，得到关节角度，角速度，角加速度，50为采样点个数


%qq=q(61,3)-q(1,3);%需要旋转的角度
%robot.plot(q);

numberTran3(arduino,q(61,1),q(61,2),q(61,3),q(61,4));%将计算结果转换为电机参数,并传入arduino

%画图
% 生成示例数据（可以替换为实际数据）
num_samples = 61; % 假设采样点数为100


% 数据转换
q(:, 1) = q(:, 1)*100; % 第一列数据单位从米变为厘米
q(:, 2:4) = q(:, 2:4) * (180 / pi); % 第二至四列数据从弧度变为度

% 创建图表
figure;
hold on;

% 绘制第一列数据，使用左纵坐标轴
yyaxis left;
plot(q(:, 1), 'LineWidth', 1.5);
ylabel('Hight (cm)');

% 绘制第二至四列数据，使用右纵坐标轴
yyaxis right;
plot(q(:, 2), 'LineWidth', 3);
plot(q(:, 3), 'LineWidth', 3);
plot(q(:, 4), 'LineWidth', 3);
ylabel('Angle (degrees)');%Angle (degrees)

% 设置横坐标轴
xlabel('Sample Points');

% 设置图表属性
title('Virtual Model Joint Data');
legend('Hight (cm)', 'J2 (degrees)','J3 (degrees)','J4 (degrees)', 'Location', 'Best');
set(gca, 'FontName', 'Times New Roman');

% 调整图表布局
grid on;
hold off;

%save('虚控实1.mat', 'q');
