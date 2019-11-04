function perturb_plot()
% this function give the plot of sqrt(phi_p/phi_w) vs q & alpha and 
% sqrt(phi_p/phi_w) vs q & delta Br

%% read data
close all;
md=importdata('../orbit_results/dbrad.plt',' ',1);
md1=importdata('../orbit_results/aharmonics.plt',' ',3);
md2=importdata('../orbit_results/profiles.plt',' ',2);
prof=md2.data;
dbrad=md.data;
aharmonics=md1.data;
xi=linspace(0,1,100);
pol=prof(:,1);
q=prof(:,3);
poln=pol/max(pol);
qi=interp1(sqrt(poln),q,xi,'PCHIP');
[m,n]=size(aharmonics);
x=aharmonics(:,1);
%% sqrt(phi_p/phi_w) vs q & alpha
figure
% vectors after position means the proportion of picture to the screen
subplot('Position',[0.075 0.15 0.3 0.8])
set(gcf,'Units','centimeters','Position',[6 6 20 15]);
for j=1:n-1
a=aharmonics(:,j+1);
ai=interp1(sqrt(x),a,xi,'PCHIP');
[ax,h1,h2]=plotyy(xi,ai,xi,qi,'plot');
hold on
end
hold off
xlabel('$\sqrt{\psi _p/\psi _w}$','FontSize',16,'interpreter','latex');
set(get(ax(1),'ylabel'),'string','$\alpha$','FontSize',16,'interpreter','latex');
set(get(ax(2),'ylabel'),'string','$q$','FontSize',16,'interpreter','latex');
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

%% sqrt(phi_p/phi_w) vs q & delta Br
subplot('Position',[0.6 0.15 0.3 0.8])
for j=1:n-1
db=dbrad(:,j+1);
dbi=interp1(sqrt(x),db,xi,'PCHIP');
[ax,h1,h2]=plotyy(xi,dbi,xi,qi);
hold on
end
hold off
xlabel('$\sqrt{\psi _p/\psi _w}$','FontSize',16,'interpreter','latex');
set(get(ax(1),'ylabel'),'string','$\delta B_{r}$','FontSize',16,'interpreter','latex');
set(get(ax(2),'ylabel'),'string','$q$','FontSize',16,'interpreter','latex');
set(ax,'FontSize',16)
set(h1,'LineWidth',3)
set(h2,'LineWidth',3)
set(ax,'LineWidth',3)
grid on
saveas(gcf,'../pictures/dbrad.fig')
saveas(gcf,'../pictures/dbrad.png')