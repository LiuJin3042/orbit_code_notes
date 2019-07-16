function FILD()
% this function gives the plot of N_loss vs t, N_loss vs ptch, N_loss in
% phase space.
% output folder is 'pictures'

%% read files and data cleaning
% creat a new dir to put pictures
% ignore the warning of exsisting files
warning off
close all
mkdir('../pictures')
% creat a file to brief the output
filename='../lost_FILD.txt';
fid=fopen(filename,'w');
m=importdata('../orbit_results/mupplane.plt',' ',1);
m1=importdata('../orbit_results/lost.plt',' ',6);% 6 rows for the column header.
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
lost=m1.data;
thet=lost(:,2);
zet=lost(:,3);
x=lost(:,4);
z=lost(:,5);
en=lost(:,6);
ptch=lost(:,7);
P=lost(:,8);
t=lost(:,13);
n=length(thet);
k=0;
t0=-pi/30; % FILD position midplane [-thet,thet], -thet=t0, thet=tm
tm=pi/30;  % or FILD position  [t0, tm]
phi0=-pi/6.0;  % FILD position midplane [phi0,phi1]
phi1=pi/6.0;   % or FILD position [-phi, phi], -phi=phi0, phi=phi1
ptch0=0.5;
en0=64;
time0=2.0;
pt = 4; 
% pt ==1 for [phi_min, phi_max] and [-thet, thet]
% pt ==2 for [phi_min, phi_max] and [thet_min, thet_max]
% pt ==3 for [-phi, phi] and [thet_min, thet_max]
% pt ==4 for [-phi, phi] and [-thet, thet]
if pt==1
   for j=1:n
    zz=rem(thet(j),2*pi);
    zz1=rem(zet(j),2*pi);
   
    %if (abs(zz1-2*pi)<=phi0 || abs(zz1)<=phi0) && (abs(zz)<=t0 || abs(zz-2*pi)<=t0) %...
   %&& t(j)>=time0 %ptch(j)>ptch0 %&& t(j)>=time0 %&& en(j)<en0 
   % obtain index of particles in FILD
   if ((zz1>=phi0 && zz1<=phi1) || (zz1>= (phi0-2*pi) && zz1<= (phi1-2*pi)))
      if (zz<=tm && zz>=t0)  || (zz>=2*pi+t0 && zz<=2*pi) ||(zz>=-2*pi && zz<=tm-2*pi) 
        k=k+1;
        t1(k)=t(j);
        lost1(k,:)=lost(j,:);
      end
   end
   end
end
if pt ==2
   for j=1:n
    zz=rem(thet(j),2*pi);
    zz1=rem(zet(j),2*pi);
   
    %if (abs(zz1-2*pi)<=phi0 || abs(zz1)<=phi0) && (abs(zz)<=t0 || abs(zz-2*pi)<=t0) %...
   %&& t(j)>=time0 %ptch(j)>ptch0 %&& t(j)>=time0 %&& en(j)<en0 
   % obtain index of particles in FILD
   if ((zz1>=phi0 && zz1<=phi1) || (zz1>= (phi0-2*pi) && zz1<= (phi1-2*pi)))
      if (zz<=tm && zz>=t0)  || (zz>=(t0-2*pi) && zz<=(tm-2*pi)) 
        k=k+1;
        t1(k)=t(j);
        lost1(k,:)=lost(j,:);
      end
   end
   end
end
if pt ==3
  for j=1:n
    zz=rem(thet(j),2*pi);
    zz1=rem(zet(j),2*pi);
   
    %if (abs(zz1-2*pi)<=phi0 || abs(zz1)<=phi0) && (abs(zz)<=t0 || abs(zz-2*pi)<=t0) %...
   %&& t(j)>=time0 %ptch(j)>ptch0 %&& t(j)>=time0 %&& en(j)<en0 
   % obtain index of particles in FILD
   if (zz1<=phi1 && zz>=phi0)  || (zz1>=(2*pi+phi0) && zz1<=2*pi) ||(zz1>=-2*pi && zz1<=(phi1-2*pi)) 
      if (zz<=tm && zz>=t0)  || (zz>=(t0-2*pi) && zz<=(tm-2*pi)) 
        k=k+1;
        t1(k)=t(j);
        lost1(k,:)=lost(j,:);
      end
   end
  end
end
if pt==4
  for j=1:n
    zz=rem(thet(j),2*pi);
    zz1=rem(zet(j),2*pi);
   
    %if (abs(zz1-2*pi)<=phi0 || abs(zz1)<=phi0) && (abs(zz)<=t0 || abs(zz-2*pi)<=t0) %...
   %&& t(j)>=time0 %ptch(j)>ptch0 %&& t(j)>=time0 %&& en(j)<en0 
   % obtain index of particles in FILD
   if (zz1<=phi1 && zz>=phi0)  || (zz1>=(2*pi+phi0) && zz1<=2*pi) ||(zz1>=-2*pi && zz1<=(phi1-2*pi)) 
      if (zz<=tm && zz>=t0)  || (zz>=(2*pi+t0) && zz<=2*pi) ||(zz>=-2*pi && zz<=(tm-2*pi)) 
        k=k+1;
        t1(k)=t(j);
        lost1(k,:)=lost(j,:);
      end
   end
  end
end
en1=lost1(:,6);
ptch1=lost1(:,7);
P1=lost1(:,8);
P10=lost1(:,9);
muB1=lost1(:,10);

%% N_loss vs t
figure(1)
num=length(en);
nt=50;
x1=linspace(min(t),max(t),nt);
h=x1(2)-x1(1);
f=zeros(1,nt);
for i=1:num
    for j=1:nt-1
        if t(i)>=x1(j) && t(i) <=x1(j+1)
            tt=t(i)-x1(j);
            f(j)=f(j)+(h-tt)/h;
            f(j+1)=f(j+1)+tt/h;
        end
    end
end
subplot(2,1,1)
%hist(t,x1)
xi=linspace(min(t),max(t),1000);
yi=interp1(x1,f,xi,'cubic');
plot(xi,yi,'LineWidth',3)
xlabel('$ t/msec$','FontSize',16,'interpreter','latex')
ylabel('$ N_{loss}$','FontSize',16,'interpreter','latex')
set(gca,'FontSize',16,'LineWidth',3)

% N_FILD vs t

subplot(2,1,2)
num1=length(en1);
nt=50;
x1=linspace(min(t1),max(t1),nt);
h1=x1(2)-x1(1);
g=zeros(1,nt);
for i=1:num1
    for j=1:nt-1
        if t1(i)>=x1(j) && t1(i) <=x1(j+1)
            tt=t1(i)-x1(j);
            g(j)=g(j)+(h-tt)/h;
            g(j+1)=g(j+1)+tt/h;
        end
    end
end

xi=linspace(min(t1),max(t1),1000);
yi=interp1(x1,g,xi,'cubic');
plot(xi,yi,'LineWidth',3);
xlabel({'$t/msec$'},'FontSize',16,'interpreter','latex')
ylabel('$ N_{loss}^{FILD}$','FontSize',16,'interpreter','latex')
set(gca,'FontSize',16,'LineWidth',3)
theta0=t0/pi*180;
theta1=tm/pi*180;
ph0=phi0/pi*180;
ph1=phi1/pi*180;
str=['\theta=' '[' num2str(theta0)  '^{o}'  ',' num2str(theta1)  '^{o}'  ']'];
str1=['\phi=' '[' num2str(ph0)  '^{o}'  ',' num2str(ph1)  '^{o}'  ']'];
title([str ', ' str1])
n2=length(lost1(:,1));
str2=['lost particle number:',num2str(n2)];
fprintf(fid,'%s\n',str2);
%%% #  lost particles-
%%% pol/pw,thet,zet,x,z,en,ptch,P,P0,mu*B,pol0,V,time(msec),type,k %%%
save('../pictures/lost1.dat','lost1','-ASCII')
saveas(gcf,'../pictures/Nloss_t.fig')
saveas(gcf,'../pictures/Nloss_t.png')

%% N_FILD in (mu*B/E, P) phase space
figure(2)
plot(P10,muB1./en1,'o',P1,muB1./en1,'rs',a,b,'m-',...
    x3,y3,'m-',x3,z3,'m-',u3,s3,'m-',v3,r3,'m-',t3,w3,'m-','LineWidth',3);
xlabel('$P_{\phi}/\psi_{w}$','FontSize',16,'interpreter','latex')
ylabel('$ \mu \cdot B _{0}/E$','FontSize',16,'interpreter','latex')
set(gca,'FontSize',16,'LineWidth',3)
legend('initial','lost')
legend boxoff
saveas(gcf,'../pictures/mu_pz.fig')
saveas(gcf,'../pictures/mu_pz.png')

%% N_FILD vs ptch
figure(3)
h=(max(ptch1)-min(ptch1))/20;
x1=min(ptch1):h:max(ptch1);
hist(ptch1,x1);
xlabel('$\lambda =\frac{v_{//}}{v}$','FontSize',16,'interpreter','latex')
ylabel('$N_{loss}$','FontSize',16,'interpreter','latex')
set(gca,'FontSize',16,'LineWidth',3)
saveas(gcf,'../pictures/Nloss_ptch.fig')
saveas(gcf,'../pictures/Nloss_ptch.png')
% N_FILD
p1=0.4;
p2=0.75;
k=0;
kk=0;
for j=1:num1
    if ptch1(j)>p1 && ptch1(j)<=p2
        k=k+1;
    else
        kk=kk+1;
    end
end

str3=['lost particle number in ptch [0.4,0.75]:',num2str(k)];
fprintf(fid,'%s\n',str3);  
str4=['lost particle number in ptch [0.75,1.0]:',num2str(kk)];
fprintf(fid,'%s\n',str4);  
fclose(fid);

