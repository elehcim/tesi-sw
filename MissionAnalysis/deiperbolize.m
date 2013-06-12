function dv = deiperbolize(x,y,vx,vy,mu,ni,ecc,a)
% Take absolute velocities and returns delta v required to make the orbit
% closed around the primary (first elliptic orbit of Galileo Mission taken 
% as a reference).
%  * (x, y) is the position of the S/C in the inertial reference system
%  * (vx, vy) is the velocity of the S/C in the inertial reference system
%  * mu is the gravitational parameter of the target planet

p = a*(1-ecc^2);
r = p/(1+ecc*cos(ni));
x_planet = r * cos(ni);
y_planet = r * sin(ni);
x_rel = x - x_planet;
y_rel = y - y_planet;

vx_planet=sqrt(mu/p)*(-sin(ni));
vy_planet=sqrt(mu/p)*(ecc+cos(ni));

vx_rel=vx-vx_planet;
vy_rel=vy-vy_planet;
v_excess=sqrt(vx_rel^2+vy_rel^2);
%dv=v_excess;

cos_ni_infty=-sign(vx_rel)*(vx/(vx_rel^2+vy_rel^2)); %TODO
e_hyp = -1/cos_ni_infty;
a_hyp=-mu/v_excess.^2;
rp_hyp=a_hyp*(1-e_hyp);
vp_hyp=sqrt(mu*(2/rp_hyp+1/a_hyp));
% Impulse at periapse similar to Galileo Mission
% rp=a*(1-e)

R_j=71492; % km
rp = 4*R_j; % perijove (km)
ra = (sqrt(419^2+290^2)/208)*100*R_j; %apojove (km)

a_ell = rp/(1-((ra-rp)/(ra+rp)));
vp_ell=sqrt(mu*(2/rp+1/a_ell));

dv=vp_hyp-vp_ell