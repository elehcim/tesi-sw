function [X,Y]=r_syn2in(x,y,nu, a,e)
% Transform the synodic position coordinate to the inertial positions
%% Distance between primaries
r = a*(1-e^2)/(1+e*cos(nu));

%% Bring to not rotating frame (Kolomaro pag. 27, formula 2.3.1)
x_prime = x * cos(nu) - y * sin(nu);
y_prime = x * sin(nu) + x * cos(nu);

%% Dimensionalize
X = x_prime*r;
Y = y_prime*r;

