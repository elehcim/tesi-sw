function [ xdot, ydot, vxdot, vydot ] = f_gpu(t,x,y,vx,vy,mu )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
xdot=vx;
ydot=vy;
vxdot=2*vy+x-((1-mu)*(x+mu))/(((x+mu)^2+y^2)^(1.5))-(mu*(x-1+mu))/(((x-1+mu)^2+y^2)^(1.5));
vydot=-2*vx+y-(1-mu)*y/(((x+mu)^2+y^2)^(1.5))-(mu*y)/(((x-1+mu)^2+y^2)^(1.5));

end

