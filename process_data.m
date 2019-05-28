clear 
close all
clc
file_name = 'lost.plt';
m=importdata(file_name,' ',6);% 6 rows for the column header.
lost=m.data;
[~,thet,zet,x,z,en,ptch,P,P0,mu_times_B,~,V,~,~] = valuation(lost);
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
    zz=abs(abs(rem(thet,2*pi)) - pi );
    zz1=abs(abs(rem(zet,2*pi)) - pi );
    good_index = zz1 > pi - phi1 & zz > pi - tm;
    lost1 = lost(good_index,:);
end

[~,thet1,zet1,x1,z1,en1,ptch1,pz1,pz10,muB1,pol0,~,t1,~,~] = valuation(lost1);
draw_figure(t1,x1,'t1','x1',filename);
% draw_figure(x1,z1,'x1','z1','lost plt x z');
% draw_figure(pz1,muB1,'pz1','mub1',filename);
% draw_figure(pz1,en1,'E/keV','Pphi','lost plt EkeV phi');
% draw_figure(ptch1,en1,'ptch','energy','lost plt ptch en');
% draw_figure(thet1,zet1,'thet1','zet1','lost plt thet zet');

          

