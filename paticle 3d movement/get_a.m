[trajectory,t_msec,x,z,pol,theta,zeta] = read_traj('traj1.plt');
theta1 = [theta(1);theta(2)];
x1 = [x(1);x(2)];
z1 = [z(1);z(2)];
theta1 = tan(theta1);
A = [-theta1,x1.*theta1]^(-1)*z1;
R = A(1)