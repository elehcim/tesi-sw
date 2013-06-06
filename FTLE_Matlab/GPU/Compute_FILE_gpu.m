function [ ftle, dphi ] = Compute_FILE_gpu( x_T, y_T, vx_T, e_T, dx, dy, dvx, de, N, filter)
%Compute_FTLE
%   This function computes Finite Iteration Lyapunov Exponents using finite
%   difference approach.
nx=size(x_T);nx=nx(1);
ny=size(y_T);ny=ny(2);
nvx=size(vx_T);nvx=nvx(3);
ne=size(e_T);ne=ne(4);
dphi=zeros(nx,ny,nvx,ne,4,4);
ftle=zeros(nx,ny,nvx,ne);
for i=2:nx-1
	for j=2:ny-1
		for k=2:nvx-1
			for l=2:ne-1
				if filter(i,j,k,l)
					dphi(i,j,k,l,1,1)=(x_T(i+1,j,k,l)-x_T(i-1,j,k,l))/(2*dx);
					dphi(i,j,k,l,1,2)=(x_T(i,j+1,k,l)-x_T(i,j-1,k,l))/(2*dy);
					dphi(i,j,k,l,1,3)=(x_T(i,j,k+1,l)-x_T(i,j,k-1,l))/(2*dvx);
					dphi(i,j,k,l,1,4)=(x_T(i,j,k,l+1)-x_T(i,j,k,l-1))/(2*de);
					
					dphi(i,j,k,l,2,1)=(y_T(i+1,j,k,l)-y_T(i-1,j,k,l))/(2*dx);
					dphi(i,j,k,l,2,2)=(y_T(i,j+1,k,l)-y_T(i,j-1,k,l))/(2*dy);
					dphi(i,j,k,l,2,3)=(y_T(i,j,k+1,l)-y_T(i,j,k-1,l))/(2*dvx);
					dphi(i,j,k,l,2,4)=(y_T(i,j,k,l+1)-y_T(i,j,k,l-1))/(2*de);
					
					dphi(i,j,k,l,3,1)=(vx_T(i+1,j,k,l)-vx_T(i-1,j,k,l))/(2*dx);
					dphi(i,j,k,l,3,2)=(vx_T(i,j+1,k,l)-vx_T(i,j-1,k,l))/(2*dy);
					dphi(i,j,k,l,3,3)=(vx_T(i,j,k+1,l)-vx_T(i,j,k-1,l))/(2*dvx);
					dphi(i,j,k,l,3,4)=(vx_T(i,j,k,l+1)-vx_T(i,j,k,l-1))/(2*de);
					
					dphi(i,j,k,l,4,1)=(e_T(i+1,j,k,l)-e_T(i-1,j,k,l))/(2*dx);
					dphi(i,j,k,l,4,2)=(e_T(i,j+1,k,l)-e_T(i,j-1,k,l))/(2*dy);
					dphi(i,j,k,l,4,3)=(e_T(i,j,k+1,l)-e_T(i,j,k-1,l))/(2*dvx);
					dphi(i,j,k,l,4,4)=(e_T(i,j,k,l+1)-e_T(i,j,k,l-1))/(2*de);
					
					ftle(i,j,k,l)=(1/abs(N))*log(norm(squeeze(dphi(i,j,k,l,:,:))));
				else
					ftle(i,j,k,l)=0;
				end
			end
		end
	end
end

