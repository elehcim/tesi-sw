function H = EnergyH(x,y,px,py,mu)
%Energy Compute the Hamiltonian
%   
H=0.5.*((px+y).^2+(py-x).^2)-Omega(x,y,mu);
end

