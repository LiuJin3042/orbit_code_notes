function final_particle_position(distf)
%% 接收distf.plt的数据, 画出粒子最终分布的统计图
thet = distf(:,2);
thet = rem(thet,2*pi);
zet = distf(:,3);
zet = rem(zet,2*pi);
x=distf(:,4);
z=distf(:,5);
en=distf(:,6);
ptch=distf(:,7);
pz=distf(:,8);
type=distf(:,13);

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

title('final distfribution of particles')
saveas(gcf,'../pictures/distf_hist.fig')
saveas(gcf,'../pictures/distf_hist.png')