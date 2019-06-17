clear 
close all
m=importdata('profiles.plt',' ',2);
prof=m.data;
pol=prof(:,1);
x=prof(:,2);
q=prof(:,3);
xn=(x-min(x))/(max(x)-min(x));
poln=pol/max(pol);
plot(sqrt(poln)*(max(x)-min(x)),q)
grid on
xlabel('r/cm')
ylabel('q')
shot=63887;
t=2.1;
str=['#',num2str(shot),'@',num2str(t),'s'];
title(str)
saveas(gcf,'qprofile.fig')
saveas(gcf,'qprofile.png')
