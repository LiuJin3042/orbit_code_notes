function init_particle_phase(lost,distf,dist,mup)
%% init particle phase marked with different colors to identify 

%% import data(if not run as a function)
% whether the particle is confined or lost
if nargin == 0
    m = importdata('../orbit_results/lost.plt',' ',5);
    lost = m.data;
    m = importdata('../orbit_results/distf.plt',' ',5);
    distf = m.data;
    m = importdata('../orbit_results/dist.plt',' ',5);
    dist = m.data;
    m = importdata('../orbit_results/mupplane.plt',' ',1);
    mup = m.data;
end

total_particle = [lost;distf]; % 'lost' now stores the whole info of particles
%% draw boundary lines with mupplane.plt
x3=mup(:,1);
y3=mup(:,2);
z3=mup(:,3);
u3=mup(:,4);
s3=mup(:,5);
v3=mup(:,6);
r3=mup(:,7);
t3=mup(:,8);
w3=mup(:,9);
hold on
a = plot(x3,y3,'k-',x3,z3,'k-',u3,s3,'k-',v3,r3,'k-',t3,w3,'k-','LineWidth',1);
set(a,'handlevisibility','off');
xlabel('$P_{\phi}/\psi_{w}$','FontSize',16,'interpreter','latex')
ylabel('$\mu \cdot B_{0}/E$','FontSize',16,'interpreter','latex')
axis([-3.5 2 0 1.4])

%% draw lost particles
lost_type = total_particle(:,14);
particle_code = total_particle(:,13);
% 损失粒子序号
lost_particle = particle_code(find(lost_type == 2 | lost_type == 4 | lost_type == 6 | lost_type == 9));
lost_init_phase = dist(lost_particle,:); % 损失粒子初始信息
lost_en = dist(lost_particle,6); % 粒子初始能量
lost_muB = dist(lost_particle,10);
lost_pz = dist(lost_particle,8);
plot(lost_pz,lost_muB./lost_en,'r.','Markersize',0.05);

conf_particle = particle_code(find(lost_type == 1 | lost_type == 3 | lost_type == 5 | lost_type == 8));
conf_init_phase = dist(conf_particle,:);
conf_en = dist(conf_particle,6); % 粒子初始能量
conf_muB = dist(conf_particle,10);
conf_pz = dist(conf_particle,8);
plot(conf_pz,conf_muB./conf_en,'g.','Markersize',0.05);

legend('lost particles','confined particles')
hold off
legend
saveas(gcf,'../pictures/init_particle_phase.fig')
saveas(gcf,'../pictures/init_particle_phase.png')











