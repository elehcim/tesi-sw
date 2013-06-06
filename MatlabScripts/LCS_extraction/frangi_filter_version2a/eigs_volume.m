function [Lambda1,Lambda2,Lambda3,Vx,Vy,Vz] = eigs_volume(Volume,Sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(Volume,Sigma);
if(nargout>2)
	[Lambda1,Lambda2,Lambda3,Vx,Vy,Vz]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
else
	[Lambda1,Lambda2,Lambda3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
end
end

