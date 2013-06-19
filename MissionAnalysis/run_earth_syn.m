t=linspace(0,2*pi,12)
sun_grav_par = 132712439935; % km^3/s^2
a_jup = 778412027; % km
ecc_jup = 0.04839266;
R = 149600000;	% km
omega= 29.7830/R;	% rad/s

for j=1:length(t)
[x(j),y(j)]=earth_synodic(t(j), sun_grav_par, a_jup, ecc_jup, omega, R);
end