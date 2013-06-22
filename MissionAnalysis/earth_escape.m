function dv=earth_escape(x,y,vx,vy)
% assuming circular earth orbit around the sun

R_earth=6378; % km
GM_earth=3.986e5; % km^3/s^2

v_planet=29.783; % km/s

% LEO orbit altitude
h = 300; % km

% Earth orbital velocity
[theta,rho]=cart2pol(x,y);

vx_earth=-v_planet*sin(theta);
vy_earth=v_planet*cos(theta);

v_excess=sqrt((vx-vx_earth)^2+(vy-vy_earth)^2);

rp=R_earth+h;

v_c=sqrt(GM_earth/rp);

a_hyp = -GM_earth/v_excess.^2;

vp_hyp=sqrt(GM_earth * (2/rp - 1/a_hyp));

dv=vp_hyp-v_c;