function E = Energy(x,y,vx,vy,mu)
%Energy Compute the Energy
%   
E=0.5.*(vx.^2+vy.^2)-Omega(x,y,mu);
end

