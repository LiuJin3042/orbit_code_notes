
close all
m=importdata('../orbit_results/distf.plt',' ',5);% 6 rows for the column header.
distf=m.data;
pol=distf(:,1);
thet=distf(:,2);
zet=distf(:,3);
x=distf(:,4);
z=distf(:,5);
en=distf(:,6);
ptch=distf(:,7);
pz=distf(:,8);
num=length(thet);
pol=x;
ptch=z;
polmin=150;
polmax=250;
ptchmin=-100;
ptchmax=100;
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

figure(1)
[X,Y]=ndgrid(poli,ptchi);
%contour(X,Y,f);
poli1=linspace(polmin,polmax,250);
ptchi1=linspace(ptchmin,ptchmax,250);
[X1,Y1]=ndgrid(poli1,ptchi1);
Z1=griddata(X,Y,f,X1,Y1);
%subplot(1,2,1)
%pcolor(X1,Y1,Z1);
%pcolor(X,Y,f);
%colormap(hsv)
%shading flat
%pcolor(X,Y,g);
%xlabel('X/cm','FontSize',16)
%ylabel('Z/cm','FontSize',16)
%set(gca,'FontSize',16)
%subplot(1,2,2)
mesh(X1,Y1,Z1);
%pcolor(X,Y,f);
%colormap(hsv)
%shading flat
%pcolor(X,Y,g);
xlabel('X/cm','FontSize',16)
ylabel('Z/cm','FontSize',16)
set(gca,'FontSize',16)
saveas(gcf,'p_p_distf.fig')
saveas(gcf,'p_p_distf.png')