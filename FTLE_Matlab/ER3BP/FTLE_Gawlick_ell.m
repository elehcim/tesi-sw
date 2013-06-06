ft=tic;
% Compute the FTLE at different time for the ER3BP
% This script takes as arguments the true anomaly "ni" at each step,
% the eccentricity, the name of the folder where to save the output figures
% and ...

%% Initializations
nx=n;
x_0_min=-0.8;
x_0_max=-0.15;
dx=(x_0_max-x_0_min)/(nx-1);
x_0=x_0_min : dx : x_0_max;

nwx=n;
wx_0_min=-2;
wx_0_max=2;
dwx=(wx_0_max-wx_0_min)/(nwx-1);
wx_0=wx_0_min : dwx : wx_0_max;

ny=3;
dy=(dx+dwx)/2;
y_0=[Y_0-dy Y_0 Y_0+dy];

ne=3;
de=dy;
e_0=[E_0-de E_0 E_0+de];

%% Setting the options for the integrator
RelTol=1e-6;AbsTol=1e-6; % A caso
% RelTol=1e-12;AbsTol=1e-12; % From Short
% RelTol=1e-13;AbsTol=1e-22; % From JD James Mireles
% RelTol=3e-14;AbsTol=1e-16; % HIGH accuracy from Ross
options=odeset('AbsTol',AbsTol,'RelTol',RelTol);%,'Events',@distance);
%options='';
%% Integration
[ x_T, y_T, vx_T, e_T, filter, delta_e ] = Integrate_FTLE_Gawlick_ell( x_0, y_0, wx_0, e_0, F, mu, ecc, ni, options);

%% Compute filter for FTLE
filter_ftle=filter;
for i=2:(nx-1)
	for j=2:(ny-1)
		for k=2:(nwx-1)
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
[ ftle, dphi ] = Compute_FTLE_Gawlick_ell( x_T, y_T, vx_T, e_T, dx, dy, dwx, de, F, filter_ftle);

%% Plot goodness of integration
% Format specific for ni
formatSpec = '%10.5'; %FIXME
% Plot
figure('Name','Delta_e')
pcolor(squeeze(delta_e(:,2,:,2))')
colorbar
shading flat
hgsave([folder_name,'Delta_e_', num2str(n),'_',num2str(ni,3),'.fig'])%FIXME
%% Plot results
figure
pcolor(x_0, wx_0, squeeze(ftle(:,2,:,2))')
colorbar
shading flat
hgsave([folder_name,'FTLE_', num2str(n),'_',num2str(ni,3),'.fig'])%FIXME
frame_time(index)=toc(ft);
fprintf('Tempo per il frame %2i: %.2f\n',index,frame_time(index))