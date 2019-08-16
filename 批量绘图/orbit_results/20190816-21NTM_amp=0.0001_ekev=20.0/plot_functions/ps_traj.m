function ps_traj()
% this function gives the figure of mu*B/E vs P_phi/psi_w
close all
muB=3.0241e+01; % mu*B0 for a particle with traj2.plt
m=importdata('../orbit_results/mupplane.plt',' ',1);
m1=importdata('../orbit_results/traj2.plt',' ',1);% 1 rows for the column header.
m3=importdata('../orbit_results/stag.plt',' ',1);
stag=m3.data;
a=stag(:,1);
b=stag(:,2);
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
traj2=m1.data;
de=traj2(:,2);
q=traj2(:,4);
ptch=traj2(:,5);
en=traj2(:,6);
pz=traj2(:,7);

plot(pz,muB./en,'.',a,b,'m-',...
    x3,y3,'m-',x3,z3,'m-',u3,s3,'m-',v3,r3,'m-',t3,w3,'m-','LineWidth',2)
xlabel('$P_{\phi}/\psi_{w}$','FontSize',16,'interpreter','latex')
ylabel('$\mu \cdot B_{0}/E$','FontSize',16,'interpreter','latex')