function [ x_T, y_T, vy_T, e_T, filter, delta_e ] = Integrate_FTLE_Masdemont( x_0, y_0, vy_0, e_0, T, mu, options)
%Integrate
%   This function performs Runge-Kutta-Fehlberg integration for given
%   initial conditions to compute FTLE
nx=length(x_0);
ny=length(y_0);
nvy=length(vy_0);
ne=length(e_0);
vx_0=zeros(nx,ny,nvy,ne);
x_T=zeros(nx,ny,nvy,ne);
y_T=zeros(nx,ny,nvy,ne);
vx_T=zeros(nx,ny,nvy,ne);
vy_T=zeros(nx,ny,nvy,ne);
e_T=zeros(nx,ny,nvy,ne);
delta_e=zeros(nx,ny,nvy,ne);
%% Look for phisically meaningful points
filter=zeros(nx,ny,nvy,ne);  %0=meaningless point 1=meaningful point
useful=ones(nx,ny,nvy,ne);
%% Integrate only useful points
useful(1,:,:,1)=0;
useful(1,:,:,3)=0;
useful(3,:,:,1)=0;
useful(3,:,:,3)=0;

%% Integrate only meaningful
h=waitbar(0,'','Name','Integration in progress, please wait!');
for i=1:nx
    toc
	for j=1:ny
		waitbar(j/ny,h,sprintf('Computing j=%i and i=%i/%i',j,i,nx));
		parfor k=1:nvy
			for l=1:ne
				if useful(i,j,k,l)
					vx_0(i,j,k,l)=sqrt(2*Omega(x_0(i),y_0(j),mu)+2*e_0(l)-vy_0(k)^2);
					if isreal(vx_0(i,j,k,l))
						filter(i,j,k,l)=1;
						
						ci=[x_0(i), y_0(j), vx_0(i,j,k,l), vy_0(k)];
						[t,Y]=ode45(@f,[0 T], ci, options, mu);
						
						if abs(t(end)) < abs(T) % Consider also negative time
							filter(i,j,k,l)=3
						end
						
						x_T(i,j,k,l)=Y(end,1);
						y_T(i,j,k,l)=Y(end,2);
						vx_T(i,j,k,l)=Y(end,3);
						vy_T(i,j,k,l)=Y(end,4);
						e_T(i,j,k,l)=0.5*(vx_T(i,j,k,l)^2+vy_T(i,j,k,l)^2)-Omega(x_T(i,j,k,l),y_T(i,j,k,l),mu);
						
						% Compute the goodness of the integration
						delta_e(i,j,k,l)=abs(e_T(i,j,k,l)-e_0(l));
					end
				end
			end
		end
	end
end
close(h);