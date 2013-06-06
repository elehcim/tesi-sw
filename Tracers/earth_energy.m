% energia della Terra
v_abs_jup = 13.07;		% km/s
v_abs_earth = 29.78;	% km/s
L = 778547200;			% km
r_earth = 149600000;	% km
%x=r_earth/L
mu = 9.537e-4;
ecc = 0.04839;
ni = 2.20;
x = -r_earth/L
y = 0;
vx = 0;
vy = v_abs_earth/v_abs_jup
Energy = 0.5*(vx^2+vy^2)-Omega(x,y,mu)/(1+ecc*cos(ni)) % Gawlik formula
% energy of the figure: -1.35