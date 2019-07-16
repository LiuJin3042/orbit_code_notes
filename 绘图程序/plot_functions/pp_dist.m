function pp_dist()
% this function gives the plot of lambda vs psi_p/psi_w
% from the initial distribution of particles dist.f

close all
%% read data
m=importdata('../orbit_results/dist.plt',' ',5);% 6 rows for the column header.
dist=m.data;
pol=dist(:,1);
thet=dist(:,2);
zet=dist(:,3);
x=dist(:,4);
z=dist(:,5);
en=dist(:,6);
ptch=dist(:,7);
pz=dist(:,8);
num=length(thet);
polmin=0.0;
polmax=1.0;
ptchmin=0.45;
ptchmax=1.0;
n1=40;
n2=40;
poli=linspace(polmin,polmax,n1);
dpol=poli(2)-poli(1);
ptchi=linspace(ptchmin,ptchmax,n2);
dptch=ptchi(2)-ptchi(1);
f=zeros(n1,n2);
for l=1:num
    %if abs(ptch1(l)-0.6)<=0.1
    for i=1:n1-1
        if pol(l)>=poli(i) && pol(l)<=poli(i+1) 
             for j=1:n2-1
                 if ptch(l)>=ptchi(j) && ptch(l)<=ptchi(j+1)
                     if pol(l)>=poli(i) && pol(l) <=poli(i)+dpol/2 ...
                             && ptch(l)>=ptchi(j) && ptch(l)<=ptchi(j)+dptch/2
                         f(i,j)=f(i,j)+1;
                         
                     elseif pol(l)>=poli(i) && pol(l) <=poli(i)+dpol/2 ...
                             && ptch(l)>=ptchi(j)+dptch/2 && ptch(l)<=ptchi(j+1)
                         f(i,j+1)=f(i,j+1)+1;
                         
                     elseif pol(l)>=poli(i)+dpol/2 && pol(l) <=poli(i+1) ...
                             && ptch(l)>=ptchi(j) && ptch(l)<=ptchi(j)+dptch/2
                            f(i+1,j)=f(i+1,j)+1;
                            
                            
                     elseif pol(l)>=poli(i)+dpol/2 && pol(l) <=poli(i+1) ...
                             && ptch(l)>=ptchi(j)+dptch/2 && ptch(l)<=ptchi(j+1)
                            f(i+1,j+1)=f(i+1,j+1)+1;
                     end
                 end
             end
        end
    %end
    end
end

%% plotting psi_p/psi_w
figure(1)
[X,Y]=ndgrid(poli,ptchi);
%contour(X,Y,f);
poli1=linspace(polmin,polmax,250);
ptchi1=linspace(ptchmin,ptchmax,250);
[X1,Y1]=ndgrid(poli1,ptchi1);
Z1=griddata(X,Y,f,X1,Y1);
pcolor(X1,Y1,Z1);
%pcolor(X,Y,f);
%colormap(hsv)
shading flat
%pcolor(X,Y,g);
xlabel('$\psi_{p}/\psi_{w}$','FontSize',16,'interpreter','latex')
ylabel('$v_{//}/v$','FontSize',16,'interpreter','latex')
set(gca,'FontSize',16)
saveas(gcf,'../pictures/p_p_dist.fig')
saveas(gcf,'../pictures/p_p_dist.png')