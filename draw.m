num=xlsread('lost.xlsx', 'A1:C90'); %读取出该区域的数据作为表格
A=num(:,1); %从B矩阵取出第一列的所有行
B=num(:,2); 
C=num(:,3); 
xx=linspace(min(A),max(A),100); %产生min(A)到max(A)均摊的50个点，目的上拟合离散点数量上的不足
yy=linspace(min(B),max(B),100); 
[xt,yt]=meshgrid(xx,yy); %做成二维网格
zt = griddata(A,B,C,xt,yt,'nearest'); %用v4点的方式进行填充
 histogram2('XBinEdges',-1:1,'YBinEdges',-2:2,'BinCounts',[1 2 3 4; 5 6 7 8])
figure
surf(xt,yt,zt) %输出结果图形
N = hist3(X,'Ctrs',{0:10:50 2000:500:5000})
shading interp


