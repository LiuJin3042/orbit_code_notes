close all
%% �����ֽ���
m=importdata('../orbit_results/mupplane.plt',' ',1);
mup=m.data;
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
plot(x3,y3,'k-',x3,z3,'k-',u3,s3,'k-',v3,r3,'k-',t3,w3,'k-','LineWidth',2)
xlabel('$P_{\phi}/\psi_{w}$','FontSize',16,'interpreter','latex')
ylabel('$\mu \cdot B_{0}/E$','FontSize',16,'interpreter','latex')

%% ��ȡdistf.plt
distf = importdata('../orbit_results/distf.plt');
distf = distf.data; % ��ȡ���ղ�������
conf_particle_code = distf(:,13);
dist = importdata('../orbit_results/dist.plt');
dist = dist.data; % ��ȡ��ʼ����
en = dist(conf_particle_code,6); % ���ӳ�ʼ����
muB = dist(conf_particle_code,10);
pz = dist(conf_particle_code,8);
plot(pz,muB./en,'.');




