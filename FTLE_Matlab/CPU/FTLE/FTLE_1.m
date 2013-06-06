clear all
tic
%% Computation parameters
n=50;
T=2;

%% Initial conditions ranges
% Assigned initial conditions
mu=0.1;
[xl1,yl1,xl2,yl2,xl3,yl3,xl4,yl4,xl5,yl5]=Lagr(mu);
C_L1=2*Omega(xl1,yl1,mu);
E_0=-C_L1/2+0.03715;
Y_0=0;

nx=n;
x_0_min=-0.8;
x_0_max=-0.15;
dx=(x_0_max-x_0_min)/(nx-1);
x_0=x_0_min : dx : x_0_max;

nvx=n;
vx_0_min=-2;
vx_0_max=2;
dvx=(vx_0_max-vx_0_min)/(nvx-1);
vx_0=vx_0_min : dvx : vx_0_max;

ny=2;
dy=(dx+dvx)/2;
y_0=[Y_0-dy Y_0 Y_0+dy];


ne=2;
de=dy;
e_0=[E_0-de E_0 E_0+de];

%% Setting the options for the integrator
% RelTol=1e-6;AbsTol=1e-6; % A caso
RelTol=1e-12;AbsTol=1e-12; % From Short
% RelTol=1e-13;AbsTol=1e-22; % From JD James Mireles
% RelTol=3e-14;AbsTol=1e-16; % HIGH accuracy from Ross
options=odeset('AbsTol',AbsTol,'RelTol',RelTol);%'Events',@distance);
%options='';
%% Integration
[ x_T, y_T, vx_T, e_T, filter, delta_e ] = Integrate_FTLE( x_0, y_0, vx_0, e_0, T, mu, options);

%% Plot goodness of integration
figure('Name','Delta_e')
pcolor(squeeze(delta_e(:,2,:,2))')

%% Compute filter for FTLE
filter_ftle=filter;
% for i=2:(nx-1)
% 	for j=2:(ny-1)
% 		for k=2:(nvx-1)
% 			for l=2:(ne-1)
% 				if filter(i,j,k,l)==0 || filter (i,j,k,l)==3
% 					filter_ftle(i,j,k,l)=0;
% 					
% 					filter_ftle(i+1,j,k,l)=0;
% 					filter_ftle(i-1,j,k,l)=0;
% 					
% 					filter_ftle(i,j+1,k,l)=0;
% 					filter_ftle(i,j-1,k,l)=0;
% 					
% 					filter_ftle(i,j,k+1,l)=0;
% 					filter_ftle(i,j,k-1,l)=0;
% 					
% 					filter_ftle(i,j,k,l+1)=0;
% 					filter_ftle(i,j,k,l-1)=0;
% 				end
% 			end
% 			
% 		end
% 	end
% end

%% Compute FTLE
[ ftle, dphi ] = Compute_FTLE( x_T, y_T, vx_T, e_T, dx, dy, dvx, de, T, filter_ftle);

%% Plot results
figure
pcolor(x_0, vx_0, squeeze(ftle(:,2,:,2))')
shading flat
toc