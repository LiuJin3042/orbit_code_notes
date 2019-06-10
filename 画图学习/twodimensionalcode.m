function createfigure(xdata1, ydata1)
%CREATEFIGURE(xdata1, ydata1)
%  XDATA1:  histogram2 x data
%  YDATA1:  histogram2 y data

%  �� MATLAB �� 10-Jun-2019 17:10:42 �Զ�����

% ���� figure
figure1 = figure;

% ���� axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% ���� histogram2
histogram2(xdata1,ydata1,'DisplayName','B vs. A','Parent',axes1,...
    'BinMethod','auto');

view(axes1,[-10.5341055880749 28.3859085847133]);
box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
