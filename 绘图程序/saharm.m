clear;
close all;
md=importdata('dbrad.plt',' ',1);
md1=importdata('aharmonics.plt',' ',3);
md2=importdata('profiles.plt',' ',2);
prof=md2.data;
dbrad=md.data;
aharmonics=md1.data;
xi=linspace(0,1,100);
pol=prof(:,1);
q=prof(:,3);
poln=pol/max(pol);
qi=interp1(sqrt(poln),q,xi,'cubic');
[m,n]=size(aharmonics);
x=aharmonics(:,1);
figure
subplot(1,2,2)
for j=1:n-1
a=aharmonics(:,j+1);
ai=interp1(sqrt(x),a,xi,'cubic');
[ax,h1,h2]=plotyy(xi,ai,xi,qi,'plot');
hold on
end
hold off
xlabel('\surd(\psi_{p}/\psi_{w})','FontSize',16);
set(get(ax(1),'ylabel'),'string','\alpha','FontSize',16);
set(get(ax(2),'ylabel'),'string','q','FontSize',16);
set(ax,'FontSize',16)
set(h1,'LineWidth',3)
set(h2,'LineWidth',3)
set(ax,'LineWidth',3)
%ymin=0;
%ymax=max(a)*1.02;
%xlim([xmin xmax])
%ylim([ymin ymax])
grid on
%[m,n]=size(dbrad);
%x=dbrad(:,1);
%subplot(1,2,2)
%for j=1:n-1
%db=dbrad(:,j+1);
%dbi=interp1(sqrt(x),db,xi,'cubic');
%[ax,h1,h2]=plotyy(xi,dbi,xi,qi);
%hold on
%end
%hold off
%xlabel('\surd(\psi_{p}/\psi_{w})','FontSize',16);
%set(get(ax(1),'ylabel'),'string','\deltaB_{r}','FontSize',16);
%set(get(ax(2),'ylabel'),'string','q','FontSize',16);
%set(ax,'FontSize',16)
%set(h1,'LineWidth',3)
%set(h2,'LineWidth',3)
%set(ax,'LineWidth',3)
%grid on
m=importdata('wall.plt',' ',2);
wall=m.data;
xw=wall(:,1);
zw=wall(:,2);
subplot(1,2,1)
plot(xw,zw,'LineWidth',3)
hold on
m=importdata('axis.plt',' ',2);
axis=m.data;
x1=axis(:,1);
y1=axis(:,2);
x2=axis(:,3);
y2=axis(:,4);
plot(x1,y1,'LineWidth',3)
plot(x2,y2,'LineWidth',3)

%m=importdata('traj1.plt',' ',1);
%traj1=m.data;
%t=traj1(:,1);
%x=traj1(:,2);
%z=traj1(:,3);
%pol=traj1(:,4);
%thet=traj1(:,5);
%zet=traj1(:,6);
%plot(x,z,'r-')
m=importdata('surf.plt',' ',1);
surf=m.data;
x01=surf(:,1);
z01=surf(:,2);
x02=surf(:,3);
z02=surf(:,4);
x03=surf(:,5);
z03=surf(:,6);
x04=surf(:,7);
z04=surf(:,8);
x05=surf(:,9);
z05=surf(:,10);
x06=surf(:,11);
z06=surf(:,12);
x07=surf(:,13);
z07=surf(:,14);
x08=surf(:,15);
z08=surf(:,16);
plot(x01,z01,'LineWidth',3)
plot(x02,z02,'LineWidth',3)
plot(x03,z03,'LineWidth',3)
plot(x04,z04,'LineWidth',3)
plot(x05,z05,'LineWidth',3)
plot(x06,z06,'LineWidth',3)
plot(x07,z07,'LineWidth',3)
plot(x08,z08,'LineWidth',3)
load spokes.plt
xx1=spokes(:,1);
zz1=spokes(:,2);
xx2=spokes(:,3);
zz2=spokes(:,4);
xx3=spokes(:,5);
zz3=spokes(:,6);
xx4=spokes(:,7);
zz4=spokes(:,8);
xx5=spokes(:,9);
zz5=spokes(:,10);
xx6=spokes(:,11);
zz6=spokes(:,12);
xx7=spokes(:,13);
zz7=spokes(:,14);
xx8=spokes(:,15);
zz8=spokes(:,16);
xx9=spokes(:,17);
zz9=spokes(:,18);
xx10=spokes(:,19);
zz10=spokes(:,20);
xx11=spokes(:,21);
zz11=spokes(:,22);
xx12=spokes(:,23);
zz12=spokes(:,24);
xx13=spokes(:,25);
zz13=spokes(:,26);
xx14=spokes(:,27);
zz14=spokes(:,28);
xx15=spokes(:,29);
zz15=spokes(:,30);
xx16=spokes(:,31);
zz16=spokes(:,32);
plot(xx1,zz1,'LineWidth',3)
plot(xx2,zz2,'LineWidth',3)
plot(xx3,zz3,'LineWidth',3)
plot(xx4,zz4,'LineWidth',3)
plot(xx5,zz5,'LineWidth',3)
plot(xx6,zz6,'LineWidth',3)
plot(xx7,zz7,'LineWidth',3)
plot(xx8,zz8,'LineWidth',3)
plot(xx9,zz9,'LineWidth',3)
plot(xx10,zz10,'LineWidth',3)
plot(xx11,zz11,'LineWidth',3)
plot(xx12,zz12,'LineWidth',3)
plot(xx13,zz13,'LineWidth',3)
plot(xx14,zz14,'LineWidth',3)
plot(xx15,zz15,'LineWidth',3)
plot(xx16,zz16,'LineWidth',3)


hold off
xlabel('X/cm','FontSize',16)
ylabel('Z/cm','FontSize',16)
set(gca,'FontSize',16,'LineWidth',3)
grid on
saveas(gcf,'xz.fig')
saveas(gcf,'xz.png')
%saveas(gcf,'dbrad.fig')
%saveas(gcf,'dbrad.png')