clear 
close all
clc
m=importdata('lost.plt',' ',6);% 6 rows for the column header.
lost=m.data;
en=lost(:,6);
ptch=lost(:,7);
pz=lost(:,8);
num=length(en);
emin=10;
emax=70;
pmin=-1.;
pmax=0.0;
n1=20;
n2=20;
ei=linspace(emin,emax,n1);
de=ei(2)-ei(1);
pi=linspace(pmin,pmax,n2);
dp=pi(2)-pi(1);
f=zeros(n1,n2);
for l=1:num
    %if abs(ptch(l)-0.6)<=0.1
    for i=1:n1-1
        if en(l)>=ei(i) && en(l)<=ei(i+1)
            ten=en(l)-ei(i);
             for j=1:n2-1
                 if pz(l)>=pi(j) && pz(l)<=pi(j+1)
                     tpz=pz(l)-pi(j);
                     f(i,j)=f(i,j)+(de-ten)*(dp-tpz)/de/dp;
                     f(i+1,j)=f(i+1,j)+ten*(dp-tpz)/de/dp;
                     f(i+1,j+1)=f(i+1,j+1)+ten*tpz/de/dp;
                     f(i,j+1)=f(i,j+1)+(de-ten)*tpz/de/dp;
                 end
             end
        end
    end
    %end
end
figure(1)
[X,Y]=ndgrid(ei,pi);
contour(X,Y,f);
%shading interp
xlabel('E/keV')
ylabel('P_{\phi}')


p1min=-0.0;
p1max=1;
n3=20;
p1i=linspace(p1min,p1max,n3);
dp1=p1i(2)-p1i(1);
g=zeros(n1,n3);
for l=1:num
    for i=1:n1-1
        if en(l)>=ei(i) && en(l)<=ei(i+1) 
            ten=en(l)-ei(i);
             for j=1:n3-1
                 if ptch(l)>=p1i(j) && ptch(l)<=p1i(j+1)
                     tptch=ptch(l)-p1i(j);
                     g(i,j)=g(i,j)+(de-ten)*(dp1-tptch)/de/dp1;
                     g(i+1,j)=g(i+1,j)+ten*(dp-tptch)/de/dp1;
                     g(i+1,j+1)=g(i+1,j+1)+ten*tptch/de/dp1;
                     g(i,j+1)=g(i,j+1)+(de-ten)*tptch/de/dp1;
                         
                     
                 end
             end
        end
    end
end
figure(2)
[X,Y]=ndgrid(ei,p1i);
pcolor(X,Y,g);
xlabel('E/keV')
ylabel('\lambda=v_{||}/v')
%shading interp
%contour(X,Y,g);
%shading interp

            

