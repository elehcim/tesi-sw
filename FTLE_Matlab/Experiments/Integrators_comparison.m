% Point definition
mu=0.1;
x_0=-0.25;
y_0=0;
vx_0=0;
vy_0=0.5;

% Integration parameters
T=2;

% CPU integration with ODE
RelTol=1e-12;AbsTol=1e-12;
options=odeset('AbsTol',AbsTol,'RelTol',RelTol);
ci=[x_0,y_0,vx_0,vy_0];
tic
[t,Y]=ode45(@f,[0 T], ci, options, mu);
ODE45_time=toc
valore_finale_ODE45=[Y(end,1),Y(end,2),Y(end,3),Y(end,4)]

% % CPU integration with RK4
% tic
% x_RK = RK4( @f, [0,T], ci, mu );
% toc

% CPU integration with RKF45
% t0=0;
% tic
% [x,y,vx,vy] = RKF45_gpu(t0, T, x_0,y_0,vx_0,vy_0, mu );
% RKF45_time=toc
% valore_finale_RKF45=[x,y,vx,vy]
% 
% Delta=valore_finale_ODE45-valore_finale_RKF45
