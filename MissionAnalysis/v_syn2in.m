function [X_dot,Y_dot]=v_syn2in(x,y,vx,vy,nu, a,e, grav_par)
% Transform the synodic velocity to the inertial velocity
% pag. 5  Gurfil and Meltzer
% "Stationkeeping: Generalization to the Elliptic Restricted Three-Body Problem"
p = a*(1-e^2);
%{
% First version
r = p/(1+e*cos(nu));
r_prime = p*e*sin(nu)/(1+e*cos(nu))^2;
nu_dot = sqrt(grav_par) / (p^(3/2))*(1+e*cos(nu))^2;
X_dot = nu_dot* ( r_prime*x + r*vx );
Y_dot = nu_dot* ( r_prime*y + r*vy );
%}

X_dot = sqrt(grav_par/p)*(e*sin(nu)*x + (1+e*cos(nu))*vx);
Y_dot = sqrt(grav_par/p)*(e*sin(nu)*y + (1+e*cos(nu))*vy);