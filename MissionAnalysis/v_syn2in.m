function [X_dot,Y_dot]=v_syn2in(x, y, vx, vy, f, a, e, GM_sun)
% Transform the synodic velocity to the inertial velocity.
%
% The axes of the inertial reference system are fixed w.r.t. the fixed
% stars and positioned as follows:
%  * The X axis is from the primaries' barycenter toward the second primary
%  pericenter.
%  * The Y axis positive toward the direction of the motion of the smaller
%  primary when it is at its pericenter.
% 
% x, y, vx, vy are all relative to the synodic frame.

p = a*(1-e^2);

r = p/(1+e*cos(f));

f_dot = sqrt(GM_sun/(p^3))*(1+e*cos(f))^2;

X_dot = f_dot * r * (e*sin(f)/(1+e*cos(f))*(x*cos(f)-y*sin(f))+...
	(vx-y)*cos(f)-(x+vy)*sin(f));
Y_dot = f_dot * r * (e*sin(f)/(1+e*cos(f))*(x*sin(f)+y*cos(f))+...
	(vx-y)*sin(f)-(x+vy)*cos(f));