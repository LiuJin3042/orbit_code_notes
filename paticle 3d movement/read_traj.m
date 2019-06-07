function [A1,A2,A3,A4,A5,A6,A7] = read_traj(file)
%A1,A2,...A7 = dt,dEE,B,q,ptch,E,Pz
traj1 = importdata(file);
traj1 = traj1.data;
for i = 1:7
    eval(['A',num2str(i),'=traj1(:,i);'])
end