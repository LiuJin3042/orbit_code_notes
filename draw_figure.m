function count_mat = draw_figure(x_data,y_data,xlabel_name,ylabel_name,filename)

x_size=40;
x_min = floor(min(x_data));
x_max = ceil(max(x_data));
x_scale = (x_max - x_min)/2;
x_min = x_min - x_scale;
x_max = x_max + x_scale;
x_index=linspace(x_min,x_max,x_size);
x_step = x_index(2)-x_index(1);

y_size=40;
y_min = floor(min(y_data));
y_max = ceil(max(y_data));
y_scale = (y_max - y_min)/2;
y_min = y_min - y_scale;
y_max = y_max + y_scale;
y_index = linspace(y_min,y_max,y_size);
y_step = y_index(2)-y_index(1);
count_mat=zeros(y_size,x_size);

position = [round((y_data-y_min)/y_step)+1,round((x_data-x_min)/x_step)+1];
xdata_len = length(x_data);
for i=1:xdata_len
    count_mat(position(i,1),position(i,2))=count_mat(position(i,1),position(i,2))+1;
end

figure
% y_index=y_index/pi*180;
% x_index = x_index/pi*180;
[X,Y]=ndgrid(x_index,y_index);
contour(X,Y,count_mat');
ylabel(ylabel_name,'FontSize',16)
xlabel(xlabel_name,'FontSize',16)
pic_title = [filename,'_',inputname(1),'_',inputname(2),'contour'];
title(pic_title)
set(gca,'FontSize',16)
saveas(gcf,[pic_title,'_contour.fig'])
saveas(gcf,[pic_title,'_contour.png'])
figure
x_list=linspace(x_min,x_max,250);
y_list=linspace(y_min,y_max,250);
[X1,Y1]=ndgrid(x_list,y_list);
Z1=griddata(X,Y,count_mat',X1,Y1);
pcolor(X1,Y1,Z1);
%pcolor(X,Y,f);
%colormap(hsv)
shading flat
%pcolor(X,Y,g);help
ylabel(ylabel_name,'FontSize',16)
xlabel(xlabel_name,'FontSize',16)
set(gca,'FontSize',16)
pic_title = [filename,'_',inputname(1),'_',inputname(2),'surf'];
title(pic_title)
%title([str ', ' str1])
colorbar
saveas(gcf,[pic_title,'_grid.fig'])
saveas(gcf,[pic_title,'_grid.png'])
%shading interp
%contour(X,Y,g);
%shading interp

end