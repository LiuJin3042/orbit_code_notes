clear
close all
filename='lost.txt';
fid=fopen(filename,'w');
shot=63887;
time=5.7; % s
E0=64; % keV
ls=0;
if ls==1
    m=importdata('../orbit_results/lost.plt',' ',6);
else
    m=importdata('../orbit_results/dist.plt',' ',5);
end
m1=importdata('../orbit_results/mupplane.plt',' ',1);
m2=importdata('../orbit_results/wall.plt',' ',2);   % 2 rows for the column header.
m3=importdata('../orbit_results/stag.plt',' ',1);
dist=m.data;
mup=m1.data;
wall=m2.data;
stag=m3.data;
a=stag(:,1);
b=stag(:,2);
xx=dist(:,4);
zz=dist(:,5);
pz=dist(:,8);
e=dist(:,6);
mub=dist(:,10);
mube=mub./e;
k=dist(:,13);
otp=dist(:,14);
x=mup(:,1);
y=mup(:,2);
z=mup(:,3);
u=mup(:,4);
s=mup(:,5);
v=mup(:,6);
r=mup(:,7);
t=mup(:,8);
w=mup(:,9);
n=length(pz);
j=0;
j1=0;
j2=0;
j3=0;
j4=0;
j5=0;
% confined
pz1=0;
pz2=0;
pz3=0;
pz4=0;
pz5=0;
for i=1:n
    if otp(i)==3
        j=j+1;
        mube1(j)=mube(i); % ctr-pass confined
        pz1(j)=pz(i);
    elseif otp(i)==1
        j1=j1+1;
        mube2(j1)=mube(i); % co-pass confined
        pz2(j1)=pz(i);
    elseif otp(i)==7
        j2=j2+1;
        mube3(j2)=mube(i);  % stagnation
        pz3(j2)=pz(i);
    elseif otp(i)==8
        j3=j3+1;
        mube4(j3)=mube(i); % conf potato
        pz4(j3)=pz(i);
    elseif otp(i)==5
        j4=j4+1;
        mube5(j4)=mube(i); % trapped confined
        pz5(j4)=pz(i);
    end
end
n1=length(pz1);
n2=length(pz2);
n3=length(pz3);
n4=length(pz4);
n5=length(pz5);
H=[];
str_legend=[];
figure(1)
hold on
if n1>1
    str1=['ctr-pass confined number:',num2str(n1)];
    fprintf(fid,'%s\n',str1);
    disp(str1)
    
    h1=plot(pz1,mube1,'g.');
    H=[H,h1];
    str_legend=[str_legend;'ctr-pass  '];
    
end
if n2>1
    str2=['co-pass confined number:',num2str(n2)];
    fprintf(fid,'%s\n',str2);
    disp(str2)
    h2=plot(pz2,mube2,'b.');
    H=[H,h2];
    str_legend=[str_legend;'co-pass   '];
end
if n3>1
    str3=['stagnation confined number:',num2str(n3)];
    fprintf(fid,'%s\n',str3);
    disp(str3)
    h3=plot(pz3,mube3,'c.');
    H=[H,h3];
    str_legend=[str_legend;'stagnation'];
end
if n4>1
    str4=['potato confined number:',num2str(n4)];
    fprintf(fid,'%s\n',str4);
    disp(str4)
    h4=plot(pz4,mube4,'r.');
    H=[H,h4];
    str_legend=[str_legend;'potato    '];
end
if n5>1
    str5=['trapped confined number:',num2str(n5)];
    fprintf(fid,'%s\n',str5);
    disp(str5)
    h5= plot(pz5,mube5,'k.');
    H=[H,h5];
    str_legend=[str_legend;'trapped   '];
end
h6=plot(a,b,'m-',x,y,'m-',x,z,'m-',u,s,'m-',v,r,'m-',t,w,'m-','LineWidth',2);
str=['E_{0}=',num2str(E0),'keV'];
str_legend = cellstr(str_legend);
j=length(H);
if j==5
    legend([H(1),H(2),H(3),H(4),H(5)],char(str_legend(1)),char(str_legend(2)),...
        char(str_legend(3)),char(str_legend(4)),char(str_legend(5)))
end
if j==4
    legend([H(1),H(2),H(3),H(4)],char(str_legend(1)),char(str_legend(2)),...
        char(str_legend(3)),char(str_legend(4)))
end
if j==3
    legend([H(1),H(2),H(3)],char(str_legend(1)),char(str_legend(2)),...
        char(str_legend(3)))
end
if j==2
    legend([H(1),H(2)],char(str_legend(1)),char(str_legend(2)))
end
if j==1
    legend(H(1),char(str_legend(1)))
end

xlabel('P_{\zeta}/\psi_{w}')
ylabel('\muB_{0}/E')
title(['shot #',num2str(shot),' @',num2str(time),'s',' for confined particles ', 'with ',str])
hold off
saveas(gcf,'ps_confined.fig')
saveas(gcf,'ps_confined.png')
j=0;
j1=0;
j2=0;
j3=0;
% lost
pz6=0;
pz7=0;
pz8=0;
pz9=0;
for i=1:n
    if otp(i)==2  % co-pass lost
        j=j+1;
        mube6(j)=mube(i);
        pz6(j)=pz(i);
        x6(j)=xx(i);
        z6(j)=zz(i);
        k6(j)=k(i);
        
    elseif otp(i)==4 % ctr-pass lost
        j1=j1+1;
        mube7(j1)=mube(i);
        pz7(j1)=pz(i);
        x7(j1)=xx(i);
        z7(j1)=zz(i);
    elseif otp(i)==6 % trapped lost
        j2=j2+1;
        mube8(j2)=mube(i);
        pz8(j2)=pz(i);
        x8(j2)=xx(i);
        z8(j2)=zz(i);
    elseif otp(i)==9 % lost potato
        j3=j3+1;
        mube9(j3)=mube(i);
        pz9(j3)=pz(i);
        x9(j3)=xx(i);
        z9(j3)=zz(i);
    end
end
n0=length(pz);
n1=length(pz6);
n2=length(pz7);
n3=length(pz8);
n4=length(pz9);
str=['shot #',num2str(shot),' @',num2str(time),'s',' for lost particles'];
str0=['total number:',num2str(n0)];
H1=[];
str_legend1=[];
figure(2)
hold on
if n1>1
    str1=['co-pass lost number:',num2str(n1)];
    fprintf(fid,'%s\n',str1);
    disp(str1)
    h7=plot(pz6,mube6,'g.');
    H1=[H1,h7];
    str_legend1=[str_legend1;'co-pass   '];
    
end
if n2>1
    str2=['ctr-pass lost number:',num2str(n2)];
    fprintf(fid,'%s\n',str2);
    disp(str2)
    h8=plot(pz7,mube7,'b.');
    H1=[H1,h8];
    str_legend1=[str_legend1;'ctr-pass  '];
end
if n3>1
    str3=['trapped lost number:',num2str(n3)];
    fprintf(fid,'%s\n',str3);
    disp(str3)
    h9=plot(pz8,mube8,'c.');
    H1=[H1,h9];
    str_legend1=[str_legend1;'trapped   '];
end
if n4>1
    str4=['lost potato number:',num2str(n4)];
    fprintf(fid,'%s\n',str4);
    disp(str4)
    h10=plot(pz9,mube9,'r.');
    H1=[H1,h10];
    str_legend1=[str_legend1;'potato    '];
end
h11=plot(a,b,'m-',x,y,'m-',x,z,'m-',u,s,'m-',v,r,'m-',t,w,'m-','LineWidth',2);
xlabel('P_{\zeta}/\psi_{w}')
ylabel('\muB_{0}/E')
str=['E_{0}=',num2str(E0),'keV'];
title(['shot #',num2str(shot),' @',num2str(time),'s',' for lost particles ', 'with ',str])

str_legend1 = cellstr(str_legend1);
j=length(H1);
if j==5
    legend([H1(1),H1(2),H1(3),H1(4),H1(5)],char(str_legend1(1)),char(str_legend1(2)),...
        char(str_legend1(3)),char(str_legend1(4)),char(str_legend1(5)))
end
if j==4
    legend([H1(1),H1(2),H1(3),H1(4)],char(str_legend1(1)),char(str_legend1(2)),...
        char(str_legend1(3)),char(str_legend1(4)))
end
if j==3
    legend([H1(1),H1(2),H1(3)],char(str_legend1(1)),char(str_legend1(2)),...
        char(str_legend1(3)))
end
if j==2
    legend([H1(1),H1(2)],char(str_legend1(1)),char(str_legend1(2)))
end
if j==1
    legend(H1(1),char(str_legend1(1)))
end
hold off
fprintf(fid,'%s\n',str);
fprintf(fid,'%s\n',str0);
disp(str0)
H2=[];
str_legend2=[];
saveas(gcf,'ps_lost.fig')
saveas(gcf,'ps_lost.png')

figure(3)
hold on
xw=wall(:,1);
zw=wall(:,2);
plot(xw,zw,'-k','LineWidth',2)

if n1>1
    str1=['co-pass lost number:',num2str(n1)];
    fprintf(fid,'%s\n',str1);
    disp(str1)
    h12=plot(x6,z6,'go');
    H2=[H2,h12];
    str_legend2=[str_legend2;'co-pass   '];
    
end
if n2>1
    str2=['ctr-pass lost number:',num2str(n2)];
    fprintf(fid,'%s\n',str2);
    disp(str2)
    h13=plot(x7,z7,'bo');
    H2=[H2,h13];
    str_legend2=[str_legend2;'ctr-pass  '];
end
if n3>1
    str3=['trapped lost number:',num2str(n3)];
    fprintf(fid,'%s\n',str3);
    disp(str3)
    h14=plot(x8,z8,'co');
    H2=[H2,h14];
    str_legend2=[str_legend2;'trapped   '];
end
if n4>1
    str4=['lost potato number:',num2str(n4)];
    fprintf(fid,'%s\n',str4);
    disp(str4)
    h15=plot(x9,z9,'ro');
    H2=[H2,h15];
    str_legend2=[str_legend2;'potato    '];
end
str_legend2 = cellstr(str_legend2);
j=length(H2);

if j==5
    legend([H2(1),H2(2),H2(3),H2(4),H2(5)],char(str_legend2(1)),char(str_legend2(2)),...
        char(str_legend2(3)),char(str_legend2(4)),char(str_legend2(5)))
end

if j==4
    legend([H2(1),H2(2),H2(3),H2(4)],char(str_legend2(1)),char(str_legend2(2)),...
        char(str_legend2(3)),char(str_legend2(4)))
end

if j==3
    legend([H2(1),H2(2),H2(3)],char(str_legend2(1)),char(str_legend2(2)),...
        char(str_legend2(3)))
end

if j==2
    legend([H2(1),H2(2)],char(str_legend2(1)),char(str_legend2(2)))
end

if j==1
    legend(H2(1),char(str_legend2(1)))
end
hold off
xlabel('X/cm')
ylabel('Z/cm')
saveas(gcf,'../pictures/xz_lost_ini.fig')
saveas(gcf,'../pictures/xz_lost_ini.png')

