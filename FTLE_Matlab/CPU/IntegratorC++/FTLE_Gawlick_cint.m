clear all
tic
% Computation parameters
n=300;
T=2;

%% Initial conditions ranges

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

ny=3;
dy=(dx+dvx)/2; % FIXME assunzione sparata a caso
y_0=[Y_0-dy Y_0 Y_0+dy];

ne=3;
de=dy;
e_0=[E_0-de E_0 E_0+de];
delta_e=zeros(nx,ny,nvx,ne);

% %% Setting the options for the integrator
% RelTol=1e-6;AbsTol=1e-6; % A caso
% % RelTol=1e-12;AbsTol=1e-12; % From Short
% % RelTol=1e-13;AbsTol=1e-22; % From JD James Mireles
% % RelTol=3e-14;AbsTol=1e-16; % HIGH accuracy from Ross
% options=odeset('AbsTol',AbsTol,'RelTol',RelTol);%,'Events',@distance);
% %options='';
%% Integration
filename='ris_int.txt';
[ x_T, y_T, vx_T, e_T] = interface(filename,nx, ny, nvx, ne);
toc
% %% Compute filter for FTLE
%TODO fare una prova scommentarlo
% filter_ftle=filter;
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
% 		end
% 	end
% end

%% Energy test
for i=1:nx
	for j=1:ny
		for k=1:nvx
			for l=1:ne
% 				vy_0(i,j,k,l)=-sqrt(2*Omega(x_0(i),y_0(j),mu)+2*e_0(l)-vx_0(k)^2);
% 				if isreal(vy_0)
					delta_e(i,j,k,l)=abs(e_T(i,j,k,l)-e_0(l));
% 				end
			end
		end
	end
end
toc

%TODO Costruisci filtro per eliminare punti accanto a punti cattivi

%% Compute FTLE
[ ftle, dphi ] = Compute_FTLE( x_T, y_T, vx_T, e_T, dx, dy, dvx, de, T);


%% Plot goodness of integration
figure('Name','Delta_e')
pcolor(squeeze(delta_e(:,2,:,2))')
colorbar
shading flat

%% Plot FTLE
figure
pcolor(x_0, vx_0, squeeze(ftle(:,2,:,2))')
colorbar
shading flat
elapsed_time=toc