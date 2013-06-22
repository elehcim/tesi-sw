function [X,Y]=r_syn2in(x,y,nu, a,e)
% Transform the synodic position coordinate (x,y) to the inertial positions.
%
% The axes of the inertial reference system are fixed w.r.t. the fixed
% stars and positioned as follows:
%  * The X axis is from the primaries' barycenter toward the second primary
%  pericenter.
%  * The Y axis positive toward the direction of the motion of the smaller
%  primary when it is at its pericenter.
% 
%  * a,e are the semimajor axis and the eccentricity of the orbit of one 
% primary with respect to the other
%  * nu is the true anomaly of the smaller primary on its orbit.
%% Distance between primaries
r = a*(1-e^2)./(1+e*cos(nu));

%% Transform to fixed frame (not rotating) (Kolomaro pag. 27, formula 2.3.1)
x_prime = x .* cos(nu) - y .* sin(nu);
y_prime = x .* sin(nu) + y .* cos(nu);

%% Dimensionalize
X = x_prime.*r;
Y = y_prime.*r;

