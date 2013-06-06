clear all
tic
%% Computation parameters
n=100;
T=-10;

%% Initial conditions ranges
% Assigned initial conditions
mu=0.0121506683;
[xl1,yl1,xl2,yl2,xl3,yl3,xl4,yl4,xl5,yl5]=Lagr(mu);
C_L1=2*Omega(xl1,yl1,mu);
C=C_L1-0.035
% C=3.1638
E_0=-C/2;

X_0=0;

nx=3;
dx=1e-5;
x_0=[X_0-dx X_0 X_0+dx];

ny=n;
y_0_min=-0.8;
y_0_max=-0.15;
dy=(y_0_max-y_0_min)/(ny-1);
y_0=y_0_min : dy : y_0_max;

nvy=n;
vy_0_min=-2;
vy_0_max=2;
dvy=(vy_0_max-vy_0_min)/(nvy-1);
vy_0=vy_0_min : dvy : vy_0_max;

ne=3;
de=1e-4;
e_0=[E_0-de E_0 E_0+de];

%% Setting the options for the integrator
RelTol=1e-6;AbsTol=1e-6; % A caso
% RelTol=1e-12;AbsTol=1e-12; % From Short
% RelTol=1e-13;AbsTol=1e-22; % From JD James Mireles
% RelTol=3e-14;AbsTol=1e-16; % HIGH accuracy from Ross
options=odeset('AbsTol',AbsTol,'RelTol',RelTol);%,'Events',@distance);
% options='';
%% Integration
[ x_T, y_T, vy_T, e_T, filter, delta_e ] = Integrate_FTLE_Masdemont( x_0, y_0, vy_0, e_0, T, mu, options);

%% Plot goodness of integration
figure('Name','Delta_e')
pcolor(squeeze(delta_e(2,:,:,2))')
shading flat

%% Compute filter for FTLE
 filter_ftle=filter;
for i=2:(nx-1)
	for j=2:(ny-1)
		for k=2:(nvy-1)
			for l=2:(ne-1)
				if filter(i,j,k,l)==0 || filter (i,j,k,l)==3
					filter_ftle(i,j,k,l)=0;
					
					filter_ftle(i+1,j,k,l)=0;
					filter_ftle(i-1,j,k,l)=0;
					
					filter_ftle(i,j+1,k,l)=0;
					filter_ftle(i,j-1,k,l)=0;
					
					filter_ftle(i,j,k+1,l)=0;
					filter_ftle(i,j,k-1,l)=0;
					
					filter_ftle(i,j,k,l+1)=0;
					filter_ftle(i,j,k,l-1)=0;
				end
			end
			
		end
	end
end

%% Compute FTLE
[ ftle, dphi ] = Compute_FTLE_Masdemont( x_T, y_T, vy_T, e_T, dx, dy, dvy, de, T, filter_ftle);

%% Plot results
figure
pcolor(y_0, vy_0, squeeze(ftle(2,:,:,2))')
shading flat
elapsed_time=toc