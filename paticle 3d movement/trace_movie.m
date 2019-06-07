clc, clear, close all

%{ 
�ó����������������˶��켣�Ķ���,��֪�����˶�ʱ����,ÿһ��ʱ�̵�����
��Ҫ�ӳ������:
    plotPosition: ������Ƶ����,����plotPosVec��������ͼƬ,��ͼƬд����Ƶ����
        plotPosVec: ������λ��,����plotCurrentTrace���켣
            plotCurrentTrace: �����ӹ켣
ԭ����:
    @author: Liang Hanpu
    @bilibili ID: С�纮��
    @date: 2019/3/3
�޸�:
    @Cartman
%}

%% ��ʼ����
addpath(genpath('./plot_function'))

%% ��������
global dt x_min x_max y_min y_max z_min z_max time_end isOutVideo;
% ����ʱ��
time_end = 10;
% ʱ����
dt = 0.05;
%���Ӹ���
N = 1;
%����ʱ��
time = 0:dt:time_end;
% ����
pos = zeros(N, 3, length(time));
%% �����ͼ������.���������Һ���ģ��.ʵ����Ӧ����ȡ����
x = cos(time);
y = sin(time);
z = time;
% ͼ��߽�
x_min = -20;
x_max = -x_min;
y_min = x_min;
y_max = -y_min;
z_min = x_min;
z_max = -z_min;
% �Ƿ������Ƶͼ��
isOutVideo = true;

%% ������ʼ
for i = 1:length(time)
    pos(1,1,i) = x(i);
    pos(1,2,i) = y(i);
    pos(1,3,i) = z(i);
end

%% ��ͼ��
% ���켣ͼ��
plotPosition(pos, time)