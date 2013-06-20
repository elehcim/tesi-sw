function [x,y]=earth_synodic(f, mu, a, e, omega, R)
%
% Assumptions:
% At t=0 the position of the earth is x = 1AU, y=0
% the argument of perihelion of both planets is 0.
% * f: Anomaly of Jupiter at which I want to compute the position of the
% Earth
% * mu: Sun grav par
% * a, e: Jupiter orbital parameters
% * omega: earth orbital mean angular speed (its orbit is circular)
% * R: 1 AU

E= 2* atan(sqrt((1-e)/(1+e))*tan(f/2)); % Eccentric anomaly Proussing formula
n=sqrt(mu/a^3); %mean motion of Jupiter around the Sun
t=(E-e*sin(E))/n; % Kepler equation
ni=omega*t; % True anomaly of the earth
X=R*cos(ni); % Assuming Earth starting at X=0
Y=R*sin(ni); % Assuming Earth starting at Y=0
r=a*(1-e^2)/(1+e*cos(f));
xi=X/r;
eta=Y/r;
x=xi*cos(f)+eta*sin(f);
y=-xi*sin(f)+eta*cos(f);