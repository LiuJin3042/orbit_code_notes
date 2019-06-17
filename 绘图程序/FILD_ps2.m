clear 
close all
clc
m=importdata('lost.plt',' ',6);% 6 rows for the column header.
lost=m.data;
thet=lost(:,2);
zet=lost(:,3);
x=lost(:,4);
z=lost(:,5);
en=lost(:,6);
ptch=lost(:,7);
pz=lost(:,8);
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
% pt ==2 for [phi_min, phi_max] and [-thet, thet]
% pt ==3 for [-phi, phi] and [-thet, thet]
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
pz1=lost1(:,8);
pz10=lost1(:,9);
muB1=lost1(:,10);
num=length(en1);
emin=55;
emax=75;
pmin=-1.;
pmax=0.0;
n1=40;
n2=40;
ei=linspace(emin,emax,n1);
de=ei(2)-ei(1);
pzi=linspace(pmin,pmax,n2);
dp=pzi(2)-pzi(1);
f=zeros(n1,n2);
for l=1:num
    %if abs(ptch1(l)-0.6)<=0.1
    for i=1:n1-1
        if en1(l)>=ei(i) && en1(l)<=ei(i+1) 
             for j=1:n2-1
                 if pz1(l)>=pzi(j) && pz1(l)<=pzi(j+1)
                     if en(l)>=ei(i) && en(l) <=ei(i)+de/2 ...
                             && pz1(l)>=pzi(j) && pz1(l)<=pzi(j)+dp/2
                         f(i,j)=f(i,j)+1;
                         
                     elseif en1(l)>=ei(i) && en1(l) <=ei(i)+de/2 ...
                             && pz1(l)>=pzi(j)+dp/2 && pz1(l)<=pzi(j+1)
                         f(i,j+1)=f(i,j+1)+1;
                         
                     elseif en1(l)>=ei(i)+de/2 && en1(l) <=ei(i+1) ...
                             && pz1(l)>=pzi(j) && pz1(l)<=pzi(j)+dp/2
                            f(i+1,j)=f(i+1,j)+1;
                            
                            
                     elseif en1(l)>=ei(i)+de/2 && en1(l) <=ei(i+1) ...
                             && pz1(l)>=pzi(j)+dp/2 && pz1(l)<=pzi(j+1)
                            f(i+1,j+1)=f(i+1,j+1)+1;
                     end
                 end
             end
        end
    %end
    end
end
figure(1)
[X,Y]=ndgrid(pzi,ei);
contour(X,Y,f');
%shading interp
ylabel('E/keV','FontSize',16)
xlabel('P_{\phi}','FontSize',16)
set(gca,'FontSize',16,'LineWidth',3)
theta0=t0/pi*180;
theta1=tm/pi*180;
ph0=phi0/pi*180;
ph1=phi1/pi*180;
str=['\theta=' '[' num2str(theta0) '^{o}'  ',' num2str(theta1) '^{o}'   ']'];
str1=['\phi=' '[' num2str(ph0) '^{o}'  ',' num2str(ph1) '^{o}'  ']'];
title([str ', ' str1])
saveas(gcf,'pz_E_FILD.fig')
saveas(gcf,'pz_E_FILD.png')

p1min=-0.0;
p1max=1;
n3=20;
p1i=linspace(p1min,p1max,n3);
dp1=p1i(2)-p1i(1);
g=zeros(n1,n3);
for l=1:num
    for i=1:n1-1
        if en1(l)>=ei(i) && en1(l)<=ei(i+1) 
             for j=1:n3-1
                 if ptch1(l)>=p1i(j) && ptch1(l)<=p1i(j+1)
                     if en1(l)>=ei(i) && en1(l) <=ei(i)+de/2 ...
                             && ptch1(l)>=p1i(j) && ptch1(l)<=p1i(j)+dp1/2
                         g(i,j)=g(i,j)+1;
                         
                     elseif en1(l)>=ei(i) && en1(l) <=ei(i)+de/2 ...
                             && ptch1(l)>=p1i(j)+dp1/2 && ptch1(l)<=p1i(j+1)
                         g(i,j+1)=g(i,j+1)+1;
                         
                     elseif en1(l)>=ei(i)+de/2 && en1(l) <=ei(i+1) ...
                             && ptch1(l)>=p1i(j) && ptch1(l)<=p1i(j)+dp1/2
                            g(i+1,j)=g(i+1,j)+1;
                            
                            
                     elseif en1(l)>=ei(i)+de/2 && en1(l) <=ei(i+1) ...
                             && ptch1(l)>=p1i(j)+dp1/2 && ptch1(l)<=p1i(j+1)
                            g(i+1,j+1)=g(i+1,j+1)+1;
                     end
                 end
             end
        end
    end
end
figure(2)
ti=acos(p1i)/pi*180;
[X,Y]=ndgrid(ti,ei);
ei1=linspace(emin,emax,250);
ti1=linspace(min(ti),max(ti),250);
[X1,Y1]=ndgrid(ti1,ei1);
Z1=griddata(X,Y,g',X1,Y1);
pcolor(X1,Y1,Z1);
%pcolor(X,Y,f);
%colormap(hsv)
shading flat
%pcolor(X,Y,g);
ylabel('E/keV','FontSize',16)
xlabel('\theta^{o}=acos(v_{||}/v)','FontSize',16)
set(gca,'FontSize',16)
title([str ', ' str1])
colorbar
saveas(gcf,'angle_E_FILD.fig')
saveas(gcf,'angle_E_FILD.png')
%shading interp
%contour(X,Y,g);
%shading interp

            

