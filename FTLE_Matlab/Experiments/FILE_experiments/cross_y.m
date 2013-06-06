function [ value,isterminal,direction ] = cross_y( ~,y,~)
%cross_y event function for computing file
%   This function detect when the solution crosses the plane y=0
value=y(2);
isterminal=0;
direction=1;
end

