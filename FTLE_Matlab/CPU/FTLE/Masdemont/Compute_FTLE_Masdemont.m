function [ ftle, dphi ] = Compute_FTLE_Masdemont( x_T, y_T, vy_T, e_T, dx, dy, dvy, de, T, filter)
%Compute_FTLE
%   This function computes Finite Time Lyapunov Exponents using finite
%   difference approach.
nx=size(x_T);nx=nx(1);
ny=size(y_T);ny=ny(2);
nvy=size(vy_T);nvy=nvy(3);
ne=size(e_T);ne=ne(4);
dphi=zeros(nx,ny,nvy,ne,4,4);
ftle=zeros(nx,ny,nvy,ne);
for i=2:(nx-1)
	for j=2:(ny-1)
		for k=2:(nvy-1)
			for l=2:(ne-1)
				if filter(i,j,k,l)
					dphi(i,j,k,l,1,1)=(x_T(i+1,j,k,l)-x_T(i-1,j,k,l))/(2*dx);
					dphi(i,j,k,l,1,2)=(x_T(i,j+1,k,l)-x_T(i,j-1,k,l))/(2*dy);
					dphi(i,j,k,l,1,3)=(x_T(i,j,k+1,l)-x_T(i,j,k-1,l))/(2*dvy);
					dphi(i,j,k,l,1,4)=(x_T(i,j,k,l+1)-x_T(i,j,k,l-1))/(2*de);
					
					dphi(i,j,k,l,2,1)=(y_T(i+1,j,k,l)-y_T(i-1,j,k,l))/(2*dx);
					dphi(i,j,k,l,2,2)=(y_T(i,j+1,k,l)-y_T(i,j-1,k,l))/(2*dy);
					dphi(i,j,k,l,2,3)=(y_T(i,j,k+1,l)-y_T(i,j,k-1,l))/(2*dvy);
					dphi(i,j,k,l,2,4)=(y_T(i,j,k,l+1)-y_T(i,j,k,l-1))/(2*de);
					
					dphi(i,j,k,l,3,1)=(vy_T(i+1,j,k,l)-vy_T(i-1,j,k,l))/(2*dx);
					dphi(i,j,k,l,3,2)=(vy_T(i,j+1,k,l)-vy_T(i,j-1,k,l))/(2*dy);
					dphi(i,j,k,l,3,3)=(vy_T(i,j,k+1,l)-vy_T(i,j,k-1,l))/(2*dvy);
					dphi(i,j,k,l,3,4)=(vy_T(i,j,k,l+1)-vy_T(i,j,k,l-1))/(2*de);
					
					dphi(i,j,k,l,4,1)=(e_T(i+1,j,k,l)-e_T(i-1,j,k,l))/(2*dx);
					dphi(i,j,k,l,4,2)=(e_T(i,j+1,k,l)-e_T(i,j-1,k,l))/(2*dy);
					dphi(i,j,k,l,4,3)=(e_T(i,j,k+1,l)-e_T(i,j,k-1,l))/(2*dvy);
					dphi(i,j,k,l,4,4)=(e_T(i,j,k,l+1)-e_T(i,j,k,l-1))/(2*de);
					
					ftle(i,j,k,l)=(1/abs(T))*log(norm(squeeze(dphi(i,j,k,l,:,:))));
				else
					ftle(i,j,k,l)=0;
				end
			end
		end
	end
end