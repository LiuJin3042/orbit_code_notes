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

%% ���·��
addpath(genpath('./plot_function'))

%%ȫ�ֱ���
global dt x_min x_max y_min y_max z_min z_max isOutVideo x_label y_label z_label;
global trace duration time_start
%% ��ȡ����
[trajectory,t_msec,x,z,pol,theta,zeta] = read_traj('traj1.plt');
d1 = x;
d2 = z*1000;


%% ��������
%����ʱ��
time = trajectory;
% ʱ����
dt = time(2) - time(1);
%��ͼ��ʼʱ��(s),����ʱ��(s)
time_start = 6;
time_end = 6.05;
%��ʼ����֡��
frame_start = ceil(time_start/dt) + 1;
frame_end = ceil(time_end/dt) + 1;
%����ʱ��
duration = time_end - time_start;
%��ͼʱ��
time = time(frame_start:frame_end);
%���Ӹ���
N = 1;
%ά��
D = 2;
% N����������ʱ�̵�D��ά����������
pos = zeros(N, D, length(time));
% ����������
x_label = 'x';
y_label = 'z';
z_label = '';
%��ʷ��������,��ֹ������൲ס��Ļ.������ʾtrace��
trace = 250;
%% �����ͼ������.���������Һ���ģ��.ʵ����Ӧ����ȡ����
%����ά�ȵ�����,������
% d1 = cos(time);
% d2 = sin(time);
% d3 = time;

% ͼ��߽�
% x_min ȡd1����Сֵ�ľ���ֵ�����ֵ�ľ���ֵ�д����һ�����෴��
x_min = -max(abs(min(d1)),abs(max(d1)));
x_max = -x_min;
y_min = -max(abs(min(d2)),abs(max(d2)));
y_max = -y_min;
z_min = x_min;
z_max = -z_min;
% �Ƿ������Ƶͼ��
isOutVideo = true;

%% ������λ�ø�ֵ
for i = 1:length(time)
    for j = 1:D
        eval(['pos(1,j,i) = d',num2str(j),'(frame_start+i);'])
    end
end

%% ��ͼ��
% ���켣ͼ��
plotPosition(pos, time)
close all