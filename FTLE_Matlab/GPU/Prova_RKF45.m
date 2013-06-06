mu=0.1;
x_0=-0.25;
y_0=0;
vx_0=0;
e_0=-1.8;
vy_0=-sqrt(2*Omega(x_0,y_0,mu)+2*e_0-vx_0^2);
tic
[x_T, y_T, vx_T, vy_T]=RKF45_gpu(0,2,x_0, y_0, vx_0, vy_0, mu);
t_rkf=toc
tic
[x_T1, y_T1, vx_T1, vy_T1]=RK4_gpu(0,2,x_0, y_0, vx_0, vy_0, mu);
t_rk=toc