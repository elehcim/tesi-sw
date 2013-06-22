function [X,Y]=r_syn2in(x,y,nu, a,e)
% Transform the synodic position coordinate to the inertial positions.
% a,e are the semimajor axis and the eccentricity of the orbit of one 
% primary with respect to the other
%% Distance between primaries
r = a*(1-e^2)/(1+e*cos(nu));

%% Transform to fixed frame (not rotating) (Kolomaro pag. 27, formula 2.3.1)
x_prime = x * cos(nu) - y * sin(nu);
y_prime = x * sin(nu) + y * cos(nu);

%% Dimensionalize
X = x_prime*r;
Y = y_prime*r;

