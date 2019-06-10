num=xlsread('lost.xlsx', 'A1:C90'); %��ȡ���������������Ϊ���
A=num(:,1); %��B����ȡ����һ�е�������
B=num(:,2); 
C=num(:,3); 
xx=linspace(min(A),max(A),100); %����min(A)��max(A)��̯��50���㣬Ŀ���������ɢ�������ϵĲ���
yy=linspace(min(B),max(B),100); 
[xt,yt]=meshgrid(xx,yy); %���ɶ�ά����
zt = griddata(A,B,C,xt,yt,'nearest'); %��v4��ķ�ʽ�������
 histogram2('XBinEdges',-1:1,'YBinEdges',-2:2,'BinCounts',[1 2 3 4; 5 6 7 8])
figure
surf(xt,yt,zt) %������ͼ��
N = hist3(X,'Ctrs',{0:10:50 2000:500:5000})
shading interp


