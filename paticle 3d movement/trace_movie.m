clc, clear, close all

%{ 
该程序用来做出粒子运动轨迹的动画,已知粒子运动时间间隔,每一个时刻的坐标
主要子程序包括:
    plotPosition: 创建视频对象,调用plotPosVec函数生成图片,将图片写入视频对象
        plotPosVec: 画粒子位置,调用plotCurrentTrace画轨迹
            plotCurrentTrace: 画粒子轨迹
原作者:
    @author: Liang Hanpu
    @bilibili ID: 小风寒呐
    @date: 2019/3/3
修改:
    @Cartman
%}

%% 添加路径
addpath(genpath('./plot_function'))

%%全局变量
global dt x_min x_max y_min y_max z_min z_max isOutVideo x_label y_label z_label;
global trace duration time_start
%% 读取数据
[trajectory,t_msec,x,z,pol,theta,zeta] = read_traj('traj1.plt');
d1 = x;
d2 = z*1000;


%% 参数设置
%采样时间
time = trajectory;
% 时间间隔
dt = time(2) - time(1);
%绘图开始时间(s),结束时间(s)
time_start = 6;
time_end = 6.05;
%起始结束帧数
frame_start = ceil(time_start/dt) + 1;
frame_end = ceil(time_end/dt) + 1;
%持续时间
duration = time_end - time_start;
%绘图时间
time = time(frame_start:frame_end);
%粒子个数
N = 1;
%维度
D = 2;
% N个粒子所有时刻的D个维度坐标坐标
pos = zeros(N, D, length(time));
% 坐标轴名称
x_label = 'x';
y_label = 'z';
z_label = '';
%历史步数设置,防止输出过多挡住屏幕.设置显示trace步
trace = 250;
%% 坐标和图像设置.这里用正弦函数模拟.实际上应当读取数据
%三个维度的数据,测试用
% d1 = cos(time);
% d2 = sin(time);
% d3 = time;

% 图像边界
% x_min 取d1的最小值的绝对值和最大值的绝对值中大的哪一个的相反数
x_min = -max(abs(min(d1)),abs(max(d1)));
x_max = -x_min;
y_min = -max(abs(min(d2)),abs(max(d2)));
y_max = -y_min;
z_min = x_min;
z_max = -z_min;
% 是否输出视频图像
isOutVideo = true;

%% 给粒子位置赋值
for i = 1:length(time)
    for j = 1:D
        eval(['pos(1,j,i) = d',num2str(j),'(frame_start+i);'])
    end
end

%% 作图区
% 做轨迹图像
plotPosition(pos, time)
close all