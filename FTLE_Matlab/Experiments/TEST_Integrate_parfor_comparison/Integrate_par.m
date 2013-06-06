function [ x_T, y_T, vx_T, e_T, filter ] = Integrate_par( x_0, y_0, vx_0, e_0, T, mu, options)
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
a=zeros(nx,ny,nvx,ne);
%% Look for phisically meaningful points
filter=zeros(nx,ny,nvx,ne);  %0=meaningless point 1=meaningful point

%% Integrate only meaningful points
h=waitbar(0,'','Name','Integration in progress, please wait!');
for i=1:nx
	c1=x_0(i);
	waitbar(i/nx,h,sprintf('Computing i=%i',i));
	for j=1:ny
		c2=y_0(j);
		for l=1:ne
			c3=e_0(l);
			parfor k=1:nvx
				vy_0(i,j,k,l)=sqrt(2*Potential(c1,c2,mu)+2*c3-vx_0(k)^2);
				if isreal(vy_0(i,j,k,l))
					filter(i,j,k,l)=1;
					ci=[c1, c2, vx_0(k), vy_0(i,j,k,l)];
					[t,Y]=ode45(@f,[0 T], ci, options, mu);
					a(i,j,k,l)=length(Y);
					x_T(i,j,k,l)=Y(a(i,j,k,l),1);
					y_T(i,j,k,l)=Y(a(i,j,k,l),2);
					vx_T(i,j,k,l)=Y(a(i,j,k,l),3);
					vy_T(i,j,k,l)=Y(a(i,j,k,l),4);
					e_T(i,j,k,l)=0.5*(vx_T(i,j,k,l)^2+vy_T(i,j,k,l)^2)-Potential(x_T(i,j,k,l),y_T(i,j,k,l),mu);
				end
			end
		end
	end
end
close(h);
