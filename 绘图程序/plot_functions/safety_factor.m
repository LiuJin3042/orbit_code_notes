function safety_factor(shot,t)
% this function plots the profile of safety factor q
% need to give shot and t(s)
close all
m=importdata('../orbit_results/profiles.plt',' ',2);
prof=m.data;
pol=prof(:,1);
x=prof(:,2);
q=prof(:,3);
xn=(x-min(x))/(max(x)-min(x));
poln=pol/max(pol);
plot(sqrt(poln)*(max(x)-min(x)),q)
grid on
xlabel('$r/cm$','FontSize',16,'interpreter','latex')
ylabel('$q$','FontSize',16,'interpreter','latex')
str=['#',num2str(shot),'@',num2str(t),'s'];
title(str)
saveas(gcf,'../pictures/q_profile.fig')
saveas(gcf,'../pictures/q_profile.png')
