function poinxz(row,col,position,wall_dir,poin_dir,name)
subplot(row,col,position)
m=importdata(wall_dir,' ',2);
wall=m.data;
xw=wall(:,1);
zw=wall(:,2);
plot(xw,zw,'r')
hold on
poincare = load(poin_dir);
p=poincare(:,1);
t=poincare(:,2);
x=poincare(:,3);
z=poincare(:,4);
r=poincare(:,5);
pz=poincare(:,6);
plot(x,z,'.','MarkerSize',1.5)
axis([125,250,-60,60])
set(gcf,'Position',[10 10 2000 2000]);
xlabel('$X$','interpreter','latex')
ylabel('$Z$','interpreter','latex')
title(name,'interpreter','none');
hold off
