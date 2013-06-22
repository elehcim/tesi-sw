function dv = deiperbolize(x,y,vx,vy,GM_jup,ni,e_jup,a_jup)
% Take absolute velocities and returns delta v required to close the orbit
% around the primary (first elliptic orbit of Galileo Mission taken as a 
% reference).
%  * (x, y) is the position of the S/C in the inertial reference system
%  * (vx, vy) is the velocity of the S/C in the inertial reference system
%  * GM_jup is the gravitational parameter of the target planet
GM_sun=132712439935; % km^3/s^2

p = a_jup*(1-e_jup^2);
r_jup = p/(1+e_jup*cos(ni));
r_jup_bar= r_jup * GM_sun/(GM_sun+GM_jup);
x_planet = r_jup_bar * cos(ni);
y_planet = r_jup_bar * sin(ni);
x_rel = x - x_planet;
y_rel = y - y_planet;
% Test Jupiter SOI is 48223000 i.e. 4.822e7
% JupSOI=48223000;
% dist=sqrt(x_rel^2+y_rel^2);
% ERR=dist-JupSOI;

vx_planet=sqrt(GM_sun/p)*(-sin(ni));
vy_planet=sqrt(GM_sun/p)*(e_jup+cos(ni));

vx_rel=vx-vx_planet;
vy_rel=vy-vy_planet;

%% Compute hyperbola orbital parameters
v_excess=sqrt(vx_rel^2+vy_rel^2);
%dv=v_excess;

h_vec=cross([x_rel y_rel 0], [vx_rel vy_rel 0]);
h=h_vec(3);

a_hyp = -GM_jup/v_excess.^2;

e_hyp = sqrt(1 - h^2/(GM_jup*a_hyp));

rp_hyp=a_hyp*(1-e_hyp);
jup_radius=71492; % km !! Equatorial Radius !!

if rp_hyp < jup_radius
	dv=NaN;
	return
end
vp_hyp=sqrt(GM_jup * (2/rp_hyp - 1/a_hyp));

%% Compute target ellipse orbital parameters
rp_target = 4*jup_radius; % perijove (km)
ra_target = (sqrt(419^2+290^2)/208)*100*jup_radius; % apojove (km)

a_target = rp_target/(1-((ra_target-rp_target)/(ra_target+rp_target)));
%vp_target = sqrt(mu_jup*(2/rp_target-1/a_target));
va_target = sqrt(GM_jup*(2/ra_target-1/a_target));

%% Compute transfer ellipse orbital parameters
rp_ell_1 = rp_hyp;
ra_ell_1 = ra_target;
a_ell_1 = rp_ell_1/(1-((ra_ell_1-rp_ell_1)/(ra_ell_1+rp_ell_1)));
vp_ell_1 = sqrt(GM_jup*(2/rp_ell_1-1/a_ell_1));
va_ell_1 = sqrt(GM_jup*(2/ra_ell_1-1/a_ell_1));

%% First burn
dv_perijove = vp_ell_1 - vp_hyp;

%% Second burn
dv_apojove = va_target - va_ell_1;

%% Total
dv = abs(dv_perijove) + abs(dv_apojove);