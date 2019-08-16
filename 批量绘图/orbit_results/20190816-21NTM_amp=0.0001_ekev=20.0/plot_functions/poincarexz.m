warning off
mkdir('../pictures')
subplot(1,3,1)
m=importdata('../orbit_results/wall.plt',' ',2);
wall=m.data;
xw=wall(:,1);
zw=wall(:,2);
plot(xw,zw,'r')
hold on
load ../orbit_results/poincare.plt
p=poincare(:,1);
t=poincare(:,2);
x=poincare(:,3);
z=poincare(:,4);
r=poincare(:,5);
pz=poincare(:,6);
plot(x,z,'.','MarkerSize',1.5)
xlabel('X')
ylabel('Z')
hold off
subplot(1,3,2)
plot(pz,t,'.','MarkerSize',1.5)
ylabel('\theta')
%xlabel('\psi_{p}/\psi_{w}')
xlabel('P_{\phi}')
subplot(1,3,3)
plot(p,t,'.','MarkerSize',1.5)
ylabel('\theta')
%xlabel('\psi_{p}/\psi_{w}')
xlabel('\psi')
saveas(gcf,'../pictures/lost_hist.fig')
saveas(gcf,'../pictures/lost_hist.png')