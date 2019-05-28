num=xlsread('xztime.xlsx', 'A1:C90'); %��ȡ���������������Ϊ���
A=num(:,1); %��B����ȡ����һ�е�������
B=num(:,2); 
C=num(:,3); 
xx=linspace(min(A),max(A),100); %����min(A)��max(A)��̯��50���㣬Ŀ���������ɢ�������ϵĲ���
yy=linspace(min(B),max(B),100); 
[xt,yt]=meshgrid(xx,yy); %���ɶ�ά����
zt = griddata(A,B,C,xt,yt,'nearest'); %��v4��ķ�ʽ�������
figure
surf(xt,yt,zt) %������ͼ��
shading interp


