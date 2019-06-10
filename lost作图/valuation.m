function [pol_pw,thet,zet,x,z,en,ptch,P,P0,mu_times_B,pol0,V,k,type] = valuation(lost)
pol_pw=lost(:,1);
thet=lost(:,2);
zet=lost(:,3);
x=lost(:,4);
z=lost(:,5);
en=lost(:,6);
ptch=lost(:,7);
P=lost(:,8);
P0=lost(:,9);
mu_times_B=lost(:,10);
pol0=lost(:,11);
V=lost(:,12);
k=lost(:,13);
type = lost(:,14);
end