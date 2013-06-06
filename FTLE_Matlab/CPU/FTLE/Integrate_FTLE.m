function [ x_T, y_T, vx_T, e_T, filter, delta_e ] = Integrate_FTLE( x_0, y_0, vx_0, e_0, T, mu, options)
%Integrate
%   This function performs Runge-Kutta-Fehlberg integration for given
%   initial conditions to compute FTLE
nx=length(x_0);
ny=length(y_0);
nvx=length(vx_0);
ne=length(e_0);
vy_0=zeros(nx,ny,nvx,ne);
x_T=zeros(nx,ny,nvx,ne);
y_T=zeros(nx,ny,nvx,ne);
vx_T=zeros(nx,ny,nvx,ne);
vy_T=zeros(nx,ny,nvx,ne);
e_T=zeros(nx,ny,nvx,ne);
delta_e=zeros(nx,ny,nvx,ne);
%% Look for phisically meaningful points
filter=zeros(nx,ny,nvx,ne);  %0=meaningless point 1=meaningful point

%% Integrate only meaningful points
h=waitbar(0,'','Name','Integration in progress, please wait!');
for i=1:nx
	waitbar(i/nx,h,sprintf('Computing i=%i',i));
	for j=1:ny
		parfor k=1:nvx
			for l=1:ne
				vy_0(i,j,k,l)=-sqrt(2*Omega(x_0(i),y_0(j),mu)+2*e_0(l)-vx_0(k)^2);
				if isreal(vy_0(i,j,k,l))
					filter(i,j,k,l)=1;
					ci=[x_0(i), y_0(j), vx_0(k), vy_0(i,j,k,l)];
					[t,Y]=ode45(@f,[0 T], ci, options, mu);
% 					a(i,j,k,l)=length(Y);
% 					if abs(t(a(i,j,k,l))) < abs(T) % Consider also negative time					time
% 						filter(i,j,k,l)=3
% 					end
					x_T(i,j,k,l)=Y(end,1);
					y_T(i,j,k,l)=Y(end,2);
					vx_T(i,j,k,l)=Y(end,3);
					vy_T(i,j,k,l)=Y(end,4);
					e_T(i,j,k,l)=0.5*(vx_T(i,j,k,l)^2+vy_T(i,j,k,l)^2)-Omega(x_T(i,j,k,l),y_T(i,j,k,l),mu);
					% Check for good integration
					delta_e(i,j,k,l)=abs(e_T(i,j,k,l)-e_0(l));
				end
			end
		end
	end
end
close(h);