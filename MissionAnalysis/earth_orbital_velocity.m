v_planet=29.783; % km/s
a_jup = 778412027; % km
e_jup = 0.04839266;
p = a_jup*(1-e_jup^2);
GM_sun=132712439935; % km^3/s^2

f=1.32;
f*180/pi
r_jup = p/(1+e_jup*cos(f));
x=-0.0245;
y=0.1973;
[theta,rho]=cart2pol(x,y);
theta*180/pi

vx_earth=-v_planet*sin(theta+f)
vy_earth=v_planet*cos(theta+f)

c=sqrt(GM_sun/p)*(1+e_jup*cos(f)) % f_dot*r
% formulas back-transformations
vx=vx_earth/c *cos(f)+vy_earth/c *sin(f)-x*(e_jup*sin(f))/(1+e_jup*cos(f))-y
vy=vy_earth/c *cos(f)-vx_earth/c *sin(f)-y*(e_jup*sin(f))/(1+e_jup*cos(f))+x
norm([vx  vy])