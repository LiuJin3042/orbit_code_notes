clear
close all
A=1e-4;
m=importdata('poink.plt',' ',2);
poink=m.data;
p=poink(:,1);
t=poink(:,2);
x=poink(:,3);
z=poink(:,4);
time=poink(:,5);
e=poink(:,6);
pz=poink(:,7);
subplot(1,3,1)
plot(x,z,'.','MarkerSize',8)
xlabel('X')
ylabel('Z')
title(['A=',num2str(A)])
subplot(1,3,2)
plot(pz,t,'.','MarkerSize',8)
xlabel('P_{\phi}')
ylabel('\theta')
subplot(1,3,3)
plot(p,t,'.','MarkerSize',8)
xlabel('\psi_{p}')
ylabel('\theta')
saveas(gcf,'poink.fig')
saveas(gcf,'poink.png')