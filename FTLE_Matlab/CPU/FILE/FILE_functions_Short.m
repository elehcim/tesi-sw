clear all
tic
%% Computation parameters
n=50;
T=80;
N=5;
%% Initial conditions ranges
% Assigned initial conditions
mu=0.012151; %Earth-Moon
[xl1,yl1,xl2,yl2,xl3,yl3,xl4,yl4,xl5,yl5]=Lagr(mu);
C_L2=2*Potential(xl2,yl2,mu);
E_0=-C_L2/2;
Y_0=0;

nx=n;
x_0_min=0.1;
x_0_max=0.6;
dx=(x_0_max-x_0_min)/(nx-1);
x_0=x_0_min : dx : x_0_max;

ny=3;
y_0_min=Y_0-dx; % Assume for similar "virtual" grid
y_0_max=Y_0+dx;
dy=(y_0_max-y_0_min)/(ny-1);
y_0=y_0_min : dy : y_0_max;

nvx=n;
vx_0_min=-0.55;
vx_0_max=-0.25;
dvx=(vx_0_max-vx_0_min)/(nvx-1);
vx_0=vx_0_min : dvx : vx_0_max;

ne=3;
e_0_min=E_0-1e-1;
e_0_max=E_0+1e-1;
de=(e_0_max-e_0_min)/(ne-1);
e_0=e_0_min : de : e_0_max;

%% Setting the options for the integrator
%RelTol=1e-12;AbsTol=1e-12; % From Short
% RelTol=1e-13;AbsTol=1e-22; % From JD James Mireles
% RelTol=3e-14;AbsTol=1e-16; % HIGH accuracy from Ross
RelTol=1e-8;AbsTol=1e-8;
options=odeset('AbsTol',AbsTol,'RelTol',RelTol,'Events',@cross_y);
%options=odeset('Events',@cross_y);
%options='';
%% Integration
[ x_T, y_T, vx_T, e_T, filter ] = Integrate_FILE( x_0, y_0, vx_0, e_0, T, N, mu, options);

%% Compute FTLE
filter_ftle=filter;
for i=2:(nx-1)
	for j=2:(ny-1)
		for k=2:(nvx-1)
			for l=2:(ne-1)
				if filter(i,j,k,l)==0
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
[ ftle, dphi ] = Compute_FILE( x_T, y_T, vx_T, e_T, dx, dy, dvx, de, T, filter_ftle);

%% Plot results
figure
pcolor(x_0, vx_0, squeeze(ftle(:,2,:,2))')
shading flat
toc