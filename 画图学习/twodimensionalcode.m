function createfigure(xdata1, ydata1)
%CREATEFIGURE(xdata1, ydata1)
%  XDATA1:  histogram2 x data
%  YDATA1:  histogram2 y data

%  由 MATLAB 于 10-Jun-2019 17:10:42 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 创建 histogram2
histogram2(xdata1,ydata1,'DisplayName','B vs. A','Parent',axes1,...
    'BinMethod','auto');

view(axes1,[-10.5341055880749 28.3859085847133]);
box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
