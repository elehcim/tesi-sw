function dv=earth_escape()
R_earth=6378; % km
GM_earth=3.986e5 % km^3/s^2

% LEO orbit altitude
h = 300; % km

v_excess=

rp=R_earth+h;

v_c=sqrt(GM_earth/rp);

a_hyp = -GM_earth/v_excess.^2;

vp_hyp=sqrt(GM_earth * (2/rp - 1/a_hyp));

dv=vp_hyp-v_c;