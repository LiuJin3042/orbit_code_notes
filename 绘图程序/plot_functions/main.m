warning off
mkdir('../pictures')
%% a quick peak of lost data to decide the position of FILD
m = importdata('../orbit_results/dist.plt');
lost = m.data;
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
xlabel('thet')
ylabel('N')

subplot(2,2,2)
histogram(zet)
xlabel('zet')
ylabel('N')

subplot(2,2,3)
plot(x,z,'.')
xlabel('x')
ylabel('z')

subplot(2,2,4)
hist3([en,ptch]);
xlabel('en')
ylabel('lambda')

saveas(gcf,'../pictures/lost_hist.fig')
saveas(gcf,'../pictures/lost_hist.png')
fprintf('please input position of FILD\n')
thet0 = input('thet0 = ');
thetm = input('thetm = ');
zet0 = input('zet0 = ');
zetm = input('zetm = ');

%% q profile
input('press enter to continue')
shot = 63887; t = 2.1;
safety_factor(shot,t)
%% FILD plot
FILD(lost,thet0,thetm,zet0,zetm)
input('press enter to continue')
FILD_ps1(lost,thet0,thetm,zet0,zetm)
input('press enter to continue')
FILD_ps2(lost,thet0,thetm,zet0,zetm)
input('press enter to continue')
lost_ps1(lost)
input('press enter to continue')
pdbrad()
input('press enter to continue')
pp_dist()
input('press enter to continue')
ppxz_dist()
input('press enter to continue')
ps_traj()
input('press enter to continue')
close all