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

%% 初始条件
addpath(genpath('./plot_function'))

%% 参数设置
global dt x_min x_max y_min y_max z_min z_max time_end isOutVideo;
% 结束时间
time_end = 10;
% 时间间隔
dt = 0.05;
%粒子个数
N = 1;
%采样时间
time = 0:dt:time_end;
% 坐标
pos = zeros(N, 3, length(time));
%% 坐标和图像设置.这里用正弦函数模拟.实际上应当读取数据
x = cos(time);
y = sin(time);
z = time;
% 图像边界
x_min = -20;
x_max = -x_min;
y_min = x_min;
y_max = -y_min;
z_min = x_min;
z_max = -z_min;
% 是否输出视频图像
isOutVideo = true;

%% 迭代开始
for i = 1:length(time)
    pos(1,1,i) = x(i);
    pos(1,2,i) = y(i);
    pos(1,3,i) = z(i);
end

%% 作图区
% 做轨迹图像
plotPosition(pos, time)