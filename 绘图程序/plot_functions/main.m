warning off
mkdir('../pictures')
%% a quick peak of initial distribution of particles
input('press enter to plot dist.plt')
m = importdata('../orbit_results/dist.plt',' ',5);
dist = m.data;
init_particle_position(dist)

%% a quick peak of final distribution of particles
try
    input('press enter to plot distf.plt')
    m = importdata('../orbit_results/distf.plt',' ',5);
    distf = m.data;
    final_particle_position(distf)
catch
    distf = [];
    input('file error, does not contain expected data, press enter to continue')
end    


try
    %% a quick peak of lost data
    input('press enter to plot lost.plt')
    m = importdata('../orbit_results/lost.plt',' ',5);
    lost = m.data;
    lost_particle_position(lost)
    %% phase of lost particles
    input('press enter to plot lost phase')
    lost_phase(lost)
catch
    lost = [];
    input('file error, does not contain expected data, press enter to continue')
end

%% q profile
input('press enter to plot q profile')
safety_factor()

%% perturbation plot
input('press enter to plot perturbation')
perturb_plot()

%% initial psi-lambda and x-z
input('press enter to plot initial psi-lambda and x-z')
ppxz_dist(dist)

input('press enter to plot trajectory')
try
    m=importdata('../orbit_results/mupplane.plt',' ',1);
    mup = m.data;
    m3=importdata('../orbit_results/stag.plt',' ',1);
    stag = m3.data;
    m1=importdata('../orbit_results/traj2.plt',' ',1);% 1 rows for the column header.
    traj2 = m1.data;
    ps_traj(mup,traj2,stag)
catch
    input('file error, does not contain expected data, press enter to continue')
end       

input('press enter to plot initial phase of particles')
init_particle_phase(lost,distf,dist,mup)

input('all done, press enter to quit')
close all


%% FILD plot(this funciton has been removed)
%thet0 = input('thet0 = ');
%thetm = input('thetm = ');
%zet0 = input('zet0 = ');
%zetm = input('zetm = ');
%FILD(lost,thet0,thetm,zet0,zetm)
%input('press enter to continue')
%FILD_ps1(lost,thet0,thetm,zet0,zetm)
%input('press enter to continue')
%FILD_ps2(lost,thet0,thetm,zet0,zetm)