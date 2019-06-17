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
subplot(1,2,1)
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
[m,n]=size(dbrad);
x=dbrad(:,1);
subplot(1,2,2)
for j=1:n-1
db=dbrad(:,j+1);
dbi=interp1(sqrt(x),db,xi,'cubic');
[ax,h1,h2]=plotyy(xi,dbi,xi,qi);
hold on
end
hold off
xlabel('\surd(\psi_{p}/\psi_{w})','FontSize',16);
set(get(ax(1),'ylabel'),'string','\deltaB_{r}','FontSize',16);
set(get(ax(2),'ylabel'),'string','q','FontSize',16);
set(ax,'FontSize',16)
set(h1,'LineWidth',3)
set(h2,'LineWidth',3)
set(ax,'LineWidth',3)
grid on
saveas(gcf,'dbrad.fig')
saveas(gcf,'dbrad.png')