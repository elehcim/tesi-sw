function [ vx,vy ] = V2vxy( nu,V )

vx=-V'*sin(nu);
vy=V'*cos(nu);

end

