function lost_particle_position(lost)
%% 接收lost.plt的数据, 画出粒子最终分布的统计图
thet = lost(:,2);
thet = rem(thet,2*pi);
zet = lost(:,3);
zet = rem(zet,2*pi);
x=lost(:,4);
z=lost(:,5);
en=lost(:,6);
ptch=lost(:,7);
pz=lost(:,8);
type=lost(:,13);

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

title('lost particles distfribution')
saveas(gcf,'../pictures/lost_hist.fig')
saveas(gcf,'../pictures/lost_hist.png')