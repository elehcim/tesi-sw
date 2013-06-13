function [dv, ERR]= deiperbolize(x,y,vx,vy,mu_jup,ni,e_jup,a_jup)
% Take absolute velocities and returns delta v required to close the orbit
% around the primary (first elliptic orbit of Galileo Mission taken as a 
% reference).
%  * (x, y) is the position of the S/C in the inertial reference system
%  * (vx, vy) is the velocity of the S/C in the inertial reference system
%  * mu is the gravitational parameter of the target planet
% FIXME ERR output should disappear, now it's only for debug
mu_sun=132712439935; % km^3/s^2

p = a_jup*(1-e_jup^2);
r_jup = p/(1+e_jup*cos(ni));
r_jup_bar= r_jup * mu_sun/(mu_sun+mu_jup);
x_planet = r_jup_bar * cos(ni);
y_planet = r_jup_bar * sin(ni);
x_rel = x - x_planet;
y_rel = y - y_planet;
% Test Jupiter SOI is 48223000 i.e. 4.822e7
JupSOI=48223000;
dist=sqrt(x_rel^2+y_rel^2);
ERR=dist-JupSOI;

vx_planet=sqrt(mu_sun/p)*(-sin(ni));
vy_planet=sqrt(mu_sun/p)*(e_jup+cos(ni));

vx_rel=vx-vx_planet;
vy_rel=vy-vy_planet;
v_excess=sqrt(vx_rel^2+vy_rel^2);
%dv=v_excess;

h_vec=cross([x_rel y_rel 0], [vx_rel vy_rel 0]);
h=h_vec(3);

a_hyp = -mu_jup/v_excess.^2;

e_hyp = sqrt(1 - h^2/(mu_jup*a_hyp));
%delta_ecc=e_hyp-1

rp_hyp=a_hyp*(1-e_hyp);
jup_radius=71492; % km !! Equatorial Radius !!

if rp_hyp < jup_radius
	dv=NaN;
	return
end
vp_hyp=sqrt(mu_jup * (2/rp_hyp - 1/a_hyp));
% Impulse at periapse similar to Galileo Mission
% rp=a*(1-e)

rp_ell = 4*jup_radius; % perijove (km)
ra_ell = (sqrt(419^2+290^2)/208)*100*jup_radius; % apojove (km)

a_ell = rp_ell/(1-((ra_ell-rp_ell)/(ra_ell+rp_ell)));
vp_ell=sqrt(mu_jup*(2/rp_ell-1/a_ell));

dv_perigee=vp_ell-vp_hyp;

