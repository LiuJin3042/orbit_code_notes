function init_particle_position(dist)
%% 接收dist.plt的数据, 画出粒子初始分布的统计图
thet = dist(:,2);
thet = rem(thet,2*pi);
zet = dist(:,3);
zet = rem(zet,2*pi);
x=dist(:,4);
z=dist(:,5);
en=dist(:,6);
ptch=dist(:,7);
pz=dist(:,8);
type=dist(:,13);

subplot(2,2,1)
histogram(thet)
xlabel('$\theta$','interpreter','latex')
ylabel('$N$','interpreter','latex')

subplot(2,2,2)
histogram(zet)
xlabel('$\zet$','interpreter','latex')
ylabel('$N$','interpreter','latex')

subplot(2,2,3)
plot(x,z,'.')
xlabel('$x$','interpreter','latex')
ylabel('$z$','interpreter','latex')

subplot(2,2,4)
hist3([en,ptch]);
xlabel('$E_n$','interpreter','latex')
ylabel('$\lambda$')

title('initial distribution of particles')
saveas(gcf,'../pictures/dist_hist.fig')
saveas(gcf,'../pictures/dist_hist.png')